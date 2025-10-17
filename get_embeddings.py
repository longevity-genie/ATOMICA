from typing import Optional, List, Dict, Any, Union
from pathlib import Path
import pickle
import json
import sys

import typer
import torch
from tqdm import tqdm
from eliot import start_action, Message

from data.dataset import PDBDataset, ProtInterfaceDataset
from models.prediction_model import PredictionModel
from models.pretrain_model import DenoisePretrainModel
from models.prot_interface_model import ProteinInterfaceModel
from trainers.abs_trainer import Trainer


app = typer.Typer(help="Generate embeddings from ATOMICA models", add_completion=False)


def load_model(
    model_ckpt: Optional[Path] = None,
    model_config: Optional[Path] = None,
    model_weights: Optional[Path] = None
) -> Union[PredictionModel, ProteinInterfaceModel, DenoisePretrainModel]:
    """Load model from checkpoint or config/weights."""
    with start_action(action_type="load_model"):
        if model_ckpt:
            Message.log(message_type="loading_checkpoint", path=str(model_ckpt))
            model = torch.load(model_ckpt)
        elif model_config and model_weights:
            Message.log(message_type="loading_from_config", config=str(model_config), weights=str(model_weights))
            with open(model_config, "r") as f:
                model_config_dict = json.load(f)
            
            model_type = model_config_dict['model_type']
            if model_type in ('PredictionModel', 'DenoisePretrainModel'):
                model = PredictionModel.load_from_config_and_weights(model_config, model_weights)
            elif model_type == 'ProteinInterfaceModel':
                model = ProteinInterfaceModel.load_from_config_and_weights(model_config, model_weights)
            else:
                raise NotImplementedError(f"Model type {model_type} not implemented")
        else:
            raise ValueError("Either model_ckpt or both model_config and model_weights must be provided")
        
        return model


def process_batch(
    model: Union[PredictionModel, DenoisePretrainModel],
    items: List[Dict[str, Any]],
    dataset: Union[PDBDataset, ProtInterfaceDataset],
    device: str = "cuda"
) -> List[Dict[str, Any]]:
    """Process a single batch of items to generate embeddings."""
    outputs = []
    for item in items:
        outputs.append({"id": item["id"]})
    
    if isinstance(dataset, ProtInterfaceDataset):
        batch_items = [item["prot_data"] for item in items]
    else:
        batch_items = [item["data"] for item in items]
    
    batch = PDBDataset.collate_fn(batch_items)
    batch = Trainer.to_device(batch, device)
    return_obj = model.infer(batch)
    
    curr_block = 0
    curr_atom = 0
    for i, item in enumerate(items):
        num_blocks = len(item["data"]["B"])
        num_atoms = len(item["data"]["A"])

        outputs[i]["graph_embedding"] = return_obj.graph_repr[i].detach().cpu().numpy()
        outputs[i]["block_embedding"] = return_obj.block_repr[curr_block: curr_block + num_blocks].detach().cpu().numpy()
        outputs[i]["atom_embedding"] = return_obj.unit_repr[curr_atom: curr_atom + num_atoms].detach().cpu().numpy()
        outputs[i]["block_id"] = item["data"]["B"]
        outputs[i]["atom_id"] = item["data"]["A"]

        curr_block += num_blocks
        curr_atom += num_atoms
    
    return outputs


def process_single_item(
    model: Union[PredictionModel, DenoisePretrainModel],
    item: Dict[str, Any],
    dataset: Union[PDBDataset, ProtInterfaceDataset],
    device: str = "cuda"
) -> Optional[Dict[str, Any]]:
    """Process a single item to generate embeddings (fallback for OOM)."""
    output = {"id": item["id"]}
    data_item = item["prot_data"] if isinstance(dataset, ProtInterfaceDataset) else item["data"]
    batch = PDBDataset.collate_fn([data_item])
    batch = Trainer.to_device(batch, device)
    return_obj = model.infer(batch)
    
    output["graph_embedding"] = return_obj.graph_repr[0].detach().cpu().numpy()
    output["block_embedding"] = return_obj.block_repr.detach().cpu().numpy()
    output["atom_embedding"] = return_obj.unit_repr.detach().cpu().numpy()
    output["block_id"] = item["data"]["B"]
    output["atom_id"] = item["data"]["A"]
    
    return output


@app.command()
def main(
    data_path: Path = typer.Option(..., help="Path to the data file either in json or pickle format"),
    output_path: Path = typer.Option(..., help="Path to save the output embeddings, should be a .pkl file"),
    model_ckpt: Optional[Path] = typer.Option(None, help="Path of the model checkpoint to load"),
    model_config: Optional[Path] = typer.Option(None, help="Path of the model config to load"),
    model_weights: Optional[Path] = typer.Option(None, help="Path of the model weights to load"),
    batch_size: int = typer.Option(4, help="Batch size for processing"),
    device: str = typer.Option("cuda", help="Device to use for inference")
) -> None:
    """Generate embeddings from ATOMICA models for input data."""
    with start_action(action_type="generate_embeddings", data_path=str(data_path), output_path=str(output_path)):
        # Load model
        model = load_model(model_ckpt, model_config, model_weights)
        
        # Prepare dataset
        if isinstance(model, ProteinInterfaceModel):
            Message.log(message_type="extracting_prot_model", info="Model is ProteinInterfaceModel, extracting prot_model")
            model = model.prot_model
            dataset = ProtInterfaceDataset(data_path)
        else:
            dataset = PDBDataset(data_path)
        
        # Convert DenoisePretrainModel to PredictionModel if needed
        if isinstance(model, DenoisePretrainModel) and not isinstance(model, PredictionModel):
            model = PredictionModel.load_from_pretrained(model_ckpt)
        
        model = model.to(device)
        
        embeddings: List[Dict[str, Any]] = []
        total_batches = len(dataset) // batch_size + 1
        
        for idx in tqdm(range(0, len(dataset), batch_size), desc="Embedding data", total=total_batches):
            items = dataset.data[idx:min(idx + batch_size, len(dataset))]
            
            with start_action(action_type="process_batch", batch_idx=idx, num_items=len(items)):
                outputs = []
                
                if "CUDA out of memory" in str(torch.cuda.memory_summary()):
                    torch.cuda.empty_cache()
                
                # Try to process the batch
                try:
                    outputs = process_batch(model, items, dataset, device)
                except RuntimeError as e:
                    if "CUDA out of memory" in str(e):
                        torch.cuda.empty_cache()
                        Message.log(message_type="oom_fallback", info="CUDA out of memory, processing items one by one")
                        
                        # Fallback: process items one by one
                        for item in items:
                            try:
                                output = process_single_item(model, item, dataset, device)
                                outputs.append(output)
                            except Exception as item_error:
                                Message.log(message_type="item_error", item_id=item['id'], error=str(item_error))
                                torch.cuda.empty_cache()
                    else:
                        raise
                
                embeddings.extend(outputs)
        
        # Save embeddings
        with start_action(action_type="save_embeddings"):
            with open(output_path, "wb") as f:
                pickle.dump(embeddings, f)
            Message.log(message_type="embeddings_saved", path=str(output_path), count=len(embeddings))


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()