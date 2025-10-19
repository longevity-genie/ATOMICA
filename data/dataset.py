# Source https://github.com/THUNLP-MT/GET

import os
import pickle
import argparse
from tqdm.contrib.concurrent import process_map
from os.path import basename, splitext
from typing import List
from collections import Counter
import gzip
import orjson
import numpy as np
import torch
import biotite.structure as bs
import biotite.structure.io.pdb as pdb

from utils.logger import print_log
from .pdb_utils import Atom, VOCAB, dist_matrix_from_coords


MODALITIES = {"PP":0, "PL":1, "Pion":2, "Ppeptide":3, "PRNA":4, "PDNA":5, "RNAL":6, "CSD":7}

class Block:
    def __init__(self, symbol: str, units: List[Atom]) -> None:
        self.symbol = symbol
        self.units = units

    def __len__(self):
        return len(self.units)
    
    def __iter__(self):
        return iter(self.units)

    @property
    def coords(self):
        return np.mean([unit.get_coord() for unit in self.units], axis=0)

    def to_data(self):
        b = VOCAB.symbol_to_idx(self.symbol)
        x, a, positions = [], [], []
        for atom in self.units:
            a.append(VOCAB.atom_to_idx(atom.get_element()))
            x.append(atom.get_coord())
            positions.append(VOCAB.atom_pos_to_idx(atom.get_pos_code()))
        block_len = len(self)
        return b, a, x, positions, block_len
        

def blocks_to_data(*blocks_list: List[List[Block]]):
    B, A, X, atom_positions, block_lengths, segment_ids = [], [], [], [], [], []
    for i, blocks in enumerate(blocks_list):
        if len(blocks) == 0:
            continue
        # global node
        cur_B = [VOCAB.symbol_to_idx(VOCAB.GLB)]
        cur_A = [VOCAB.get_atom_global_idx()]
        cur_X = [None]
        cur_atom_positions = [VOCAB.get_atom_pos_global_idx()]
        cur_block_lengths = [1]
        # other nodes
        for block in blocks:
            b, a, x, positions, block_len = block.to_data()
            cur_B.append(b)
            cur_A.extend(a)
            cur_X.extend(x)
            cur_atom_positions.extend(positions)
            cur_block_lengths.append(block_len)
        # update coordinates of the global node to the center
        cur_X[0] = np.mean(cur_X[1:], axis=0).tolist()
        for x_i, x in enumerate(cur_X):
            if isinstance(x, np.ndarray):
                cur_X[x_i] = x.tolist()
        cur_segment_ids = [i for _ in cur_B]
        
        # finish these blocks
        B.extend(cur_B)
        A.extend(cur_A)
        X.extend(cur_X)
        atom_positions.extend(cur_atom_positions)
        block_lengths.extend(cur_block_lengths)
        segment_ids.extend(cur_segment_ids)

    data = {
        'X': X,   # [Natom, 2, 3]
        'B': B,             # [Nb], block (residue) type
        'A': A,             # [Natom]
        'atom_positions': atom_positions,  # [Natom]
        'block_lengths': block_lengths,  # [Nresidue]
        'segment_ids': segment_ids,      # [Nresidue]
    }

    return data


def data_to_blocks(data, fragmentation_method=None):
    if fragmentation_method:
        VOCAB.load_tokenizer(fragmentation_method)
    curr_atom_idx = 0
    list_of_blocks = []
    curr_segment_id = 0
    curr_blocks = []
    for block_idx, block in enumerate(data['B']):
        symbol = VOCAB.idx_to_symbol(block)
        if symbol == VOCAB.GLB:
            curr_atom_idx += data['block_lengths'][block_idx]
            continue
        atom_coords = data['X'][curr_atom_idx:curr_atom_idx+data['block_lengths'][block_idx]]
        atom_positions = data['atom_positions'][curr_atom_idx:curr_atom_idx+data['block_lengths'][block_idx]]
        atoms = []
        for i, atom in enumerate(data['A'][curr_atom_idx:curr_atom_idx+data['block_lengths'][block_idx]]):
            atom_name=VOCAB.idx_to_atom(atom)
            if atom_name == VOCAB.atom_global:
                continue
            element=VOCAB.idx_to_atom(atom)
            coordinate=atom_coords[i]
            pos_code=VOCAB.idx_to_atom_pos(atom_positions[i])
            atoms.append(Atom(atom_name=atom_name, element=element, coordinate=coordinate, pos_code=pos_code))
        curr_atom_idx += data['block_lengths'][block_idx]
        if data['segment_ids'][block_idx] != curr_segment_id:
            list_of_blocks.append(curr_blocks)
            curr_blocks = []
            curr_segment_id = data['segment_ids'][block_idx]
        curr_blocks.append(Block(symbol, atoms))
    list_of_blocks.append(curr_blocks)
    return list_of_blocks


def blocks_to_coords(blocks: List[Block]):
    max_n_unit = 0
    coords, masks = [], []
    for block in blocks:
        coords.append([unit.get_coord() for unit in block.units])
        max_n_unit = max(max_n_unit, len(coords[-1]))
        masks.append([1 for _ in coords[-1]])
    
    for i in range(len(coords)):
        num_pad =  max_n_unit - len(coords[i])
        coords[i] = coords[i] + [[0, 0, 0] for _ in range(num_pad)]
        masks[i] = masks[i] + [0 for _ in range(num_pad)]
    
    return np.array(coords), np.array(masks).astype('bool')  # [N, M, 3], [N, M], M == max_n_unit, in mask 0 is for padding


def df_to_blocks(df, key_residue='residue', key_insertion_code='insertion_code', key_resname='resname',
                     key_atom_name='atom_name', key_element='element', key_x='x', key_y='y', key_z='z', 
                     key_chain='chain', return_res_seq=False) -> List[Block]:
    last_res_id, last_res_symbol, last_residue = None, None, None
    blocks, units, res_seq = [], [], []
    for row in df.itertuples():  # each row is an atom (unit)
        residue = getattr(row, key_residue)
        if key_insertion_code is None:
            res_id = str(residue)
        else:
            insert_code = getattr(row, key_insertion_code)
            res_id = f'{residue}{insert_code}'.rstrip()
        if res_id != last_res_id:  # one block ended
            # if last_res_symbol == VOCAB.UNK:
            #     print('unk')
            #     print([str(a) for a in units])
            block = Block(last_res_symbol, units)
            blocks.append(block)
            res_seq.append(last_residue)
            # clear
            units = []
            last_res_id = res_id
            last_res_symbol = VOCAB.abrv_to_symbol(getattr(row, key_resname))
            last_residue = f'{getattr(row, key_chain)}_{residue}'
        atom = getattr(row, key_atom_name)
        element = getattr(row, key_element)
        if element == 'H':
            continue
        units.append(Atom(atom, [getattr(row, axis) for axis in [key_x, key_y, key_z]], element))
    blocks = blocks[1:]
    res_seq = res_seq[1:]
    blocks.append(Block(last_res_symbol, units))
    res_seq.append(last_residue)
    if return_res_seq:
        return blocks, res_seq
    return blocks


def blocks_interface(blocks1, blocks2, dist_th, return_indexes=False):
    blocks_coord, blocks_mask = blocks_to_coords(blocks1 + blocks2)
    blocks1_coord, blocks1_mask = blocks_coord[:len(blocks1)], blocks_mask[:len(blocks1)]
    blocks2_coord, blocks2_mask = blocks_coord[len(blocks1):], blocks_mask[len(blocks1):]
    dist = dist_matrix_from_coords(blocks1_coord, blocks1_mask, blocks2_coord, blocks2_mask)
    
    on_interface = dist < dist_th
    indexes1 = np.nonzero(on_interface.sum(axis=1) > 0)[0]
    indexes2 = np.nonzero(on_interface.sum(axis=0) > 0)[0]

    blocks1 = [blocks1[i] for i in indexes1]
    blocks2 = [blocks2[i] for i in indexes2]

    if return_indexes:
        return blocks1, blocks2, indexes1, indexes2
    else:
        return blocks1, blocks2


def item_to_pdb_file(item, output_pdb_file):
    atoms_list = []
    elements_list = []
    chains_list = []
    atom_names = []
    res_ids = []
    res_names = []
    hetero_flags = []

    start_atom = 0

    for block_idx, block_length in enumerate(item['data']['block_lengths']):
        if item['data']['B'][block_idx] == VOCAB.symbol_to_idx(VOCAB.GLB):
            start_atom += block_length
            continue
        for atom_idx in range(start_atom, start_atom+block_length):
            atoms_list.append(item['data']['X'][atom_idx])
            elements_list.append(VOCAB.idx_to_atom(item['data']['A'][atom_idx]))
            atom_pos = VOCAB.idx_to_atom_pos(item['data']['atom_positions'][atom_idx])
            if atom_pos == 'sm':
                atom_names.append(f"{VOCAB.idx_to_atom(item['data']['A'][atom_idx])}")
            else:
                atom_names.append(f"{VOCAB.idx_to_atom(item['data']['A'][atom_idx])}{atom_pos}")
            if block_idx in item['block_to_pdb_indexes']:
                res_ids.append(int(item['block_to_pdb_indexes'][block_idx].split("_")[1]))
                res_names.append(VOCAB.idx_to_abrv(item['data']['B'][block_idx]))
                chains_list.append(item['block_to_pdb_indexes'][block_idx].split("_")[0])
                hetero_flags.append(False)
            else:
                res_ids.append(0)
                res_names.append("UNK")
                chains_list.append("X")
                hetero_flags.append(True)
        start_atom += block_length

    coords = np.array(atoms_list)
    elements = np.array(elements_list)
    chains = np.array(chains_list)
    atom_names = np.array(atom_names)
    res_ids = np.array(res_ids)
    res_names = np.array(res_names)
    hetero_flags = np.array(hetero_flags)

    atom_array = bs.AtomArray(coords.shape[0])
    atom_array.coord = coords
    atom_array.element = elements
    atom_array.chain_id = chains
    atom_array.atom_name = atom_names
    atom_array.res_id = res_ids
    atom_array.res_name = res_names
    atom_array.hetero = hetero_flags

    pdb_file = pdb.PDBFile()
    pdb_file.set_structure(atom_array)
    with open(output_pdb_file, "w") as file:
        pdb_file.write(file)

class BlockGeoAffDataset(torch.utils.data.Dataset):

    def __init__(self, data_file, database=None, dist_th=6, n_cpu=4, suffix=''):
        '''
        data_file: path to the dataset file, can be index file in some occasions
        database: path/directory containing the complete data
        dist_th: threshold for deciding the interacting environment (minimum distance between heavy atoms of residues)
        n_cpu: number of cpus used in parallel preprocessing
        '''
        super().__init__()
        self.dist_th = dist_th
        self.data_file = os.path.abspath(data_file)
        self.database = database
        proc_file = os.path.join(
            os.path.split(data_file)[0],
            basename(splitext(data_file)[0]) + f'.{type(self).__name__}{suffix}_processed.pkl'
        )
        self.proc_file = proc_file
        need_process = True
        if os.path.exists(proc_file):
            print_log(f'Loading preprocessed data from {proc_file}...')
            with open(proc_file, 'rb') as fin:
                th, indexes, data = pickle.load(fin)
            if th == dist_th:
                self.indexes = indexes
                self.data = data
                need_process = False
        if need_process:
            print_log('Preprocessing...')
            items = self._load_data_file()
            if isinstance(items, list):
                data = process_map(self._preprocess, items, max_workers=n_cpu, chunksize=10)
            else:  # LMDB
                print('Data not list, disable parallel processing')
                from tqdm import tqdm
                data = [self._preprocess(item) for i, item in enumerate(tqdm(items))]
            data = [item for item in data if item is not None]
            self.indexes, self.data = self._post_process(items, data)
            with open(proc_file, 'wb') as fout:
                pickle.dump((dist_th, self.indexes, self.data), fout)
            print_log(f'Preprocessed data saved to {proc_file}')
    
    def _load_data_file(self):
        with open(self.data_file, 'rb') as fin:
            items = pickle.load(fin)
        return items
    
    def _post_process(self, items, processed_data):
        indexes = [ { 'id': item['id'], 'affinity': item['affinity'] } for item, d in zip(items, processed_data) if d is not None ]
        data = [d for d in processed_data if d is not None]
        return indexes, data

    def _preprocess(self, item):
        blocks1 = df_to_blocks(item['atoms_interface1'], key_atom_name='name')
        blocks2 = df_to_blocks(item['atoms_interface2'], key_atom_name='name')

        data = blocks_to_data(blocks1, blocks2)

        data['label'] = item['affinity']['neglog_aff']

        return data
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        '''
        an example of the returned data
        {
            'X': [Natom, 3],
            'B': [Nblock],
            'A': [Natom],
            'atom_positions': [Natom],
            'block_lengths': [Nblock]
            'segment_ids': [Nblock]
            'label': [1]
        }
        '''
        item = self.data[idx]
        return item
    
    @classmethod
    def filter_for_segment(cls, data, keep_segment):
        # segment_id = 0, protein
        # segment_id = 1, ligand
        for k, v in data.items():
            if type(v) is list:
                data[k] = np.array(v)

        block_mask = data['segment_ids'] == keep_segment
        block_id = np.zeros_like(data["A"]) # [Nu]
        block_id[np.cumsum(data["block_lengths"], axis=0)[:-1]] = 1
        block_id = np.cumsum(block_id, axis=0)
        atom_mask = data["segment_ids"][block_id] == keep_segment

        new_data = {}
        new_data['X'] = data['X'][atom_mask].tolist()
        new_data['B'] = data['B'][block_mask].tolist()
        new_data['A'] = data['A'][atom_mask].tolist()
        new_data['atom_positions'] = data['atom_positions'][atom_mask].tolist()
        new_data['block_lengths'] = data['block_lengths'][block_mask].tolist()
        new_data['segment_ids'] = data['segment_ids'][block_mask].tolist()
        for key in data:
            if key in ['X', 'B', 'A', 'atom_positions', 'block_lengths', 'segment_ids']:
                continue
            new_data[key] = data[key]
        return new_data

    @classmethod
    def collate_fn(cls, batch):
        keys = ['X', 'B', 'A', 'atom_positions', 'block_lengths', 'segment_ids']
        types = [torch.float, torch.long, torch.long, torch.long, torch.long, torch.long]
        res = {}
        for key, _type in zip(keys, types):
            val = []
            for item in batch:
                val.append(torch.tensor(item[key], dtype=_type))
            res[key] = torch.cat(val, dim=0)
        res['label'] = torch.tensor([item['label'] for item in batch], dtype=torch.float)
        lengths = [len(item['B']) for item in batch]
        res['lengths'] = torch.tensor(lengths, dtype=torch.long)
        return res
    

class PDBBindBenchmark(torch.utils.data.Dataset):

    def __init__(self, data_file):
        super().__init__()
        self.data = open_data_file(data_file)
        self.indexes = [ {'id': item['id'], 'label': item['affinity']['neglog_aff'] } for item in self.data ]  # to satify the requirements of inference.py

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        '''
        an example of the returned data
        {
            'X': [Natom, 3],
            'B': [Nblock],
            'A': [Natom],
            'atom_positions': [Natom],
            'block_lengths': [Natom]
            'segment_ids': [Nblock]
            'label': [1]
        }        
        '''
        item = self.data[idx]
        data = item['data']
        data['label'] = item['affinity']['neglog_aff']

        return data

    @classmethod
    def collate_fn(cls, batch):
        keys = ['X', 'B', 'A', 'block_lengths', 'segment_ids']
        types = [torch.float, torch.long, torch.long, torch.long, torch.long]
        has_block_embeddings = 'block_embeddings' in batch[0]
        has_block_embeddings_separate = 'block_embeddings0' in batch[0] and 'block_embeddings1' in batch[0]
        if has_block_embeddings:
            keys.append('block_embeddings')
            types.append(torch.float)
        elif has_block_embeddings_separate:
            keys.append('block_embeddings0')
            keys.append('block_embeddings1')
            types.append(torch.float)
            types.append(torch.float)
        res = {}
        for key, _type in zip(keys, types):
            val = []
            for item in batch:
                val.append(torch.tensor(item[key], dtype=_type))
            res[key] = torch.cat(val, dim=0)
        res['label'] = torch.tensor([item['label'] for item in batch], dtype=torch.float)
        lengths = [len(item['B']) for item in batch]
        res['lengths'] = torch.tensor(lengths, dtype=torch.long)
        return res


class PDBDataset(torch.utils.data.Dataset):

    def __init__(self, data_file: str, start_idx: int = 0, num_lines: int | None = None):
        super().__init__()
        self.data = open_data_file(data_file, start_idx=start_idx, num_lines=num_lines)
        self.indexes = [ item['id'] for item in self.data ]  # to satify the requirements of inference.py

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        '''
        an example of the returned data
        {
            'X': [Natom, 3],
            'B': [Nblock],
            'A': [Natom],
            'atom_positions': [Natom],
            'block_lengths': [Natom]
            'segment_ids': [Nblock]
        }        
        '''
        item = self.data[idx]
        data = item['data']

        return data

    @classmethod
    def collate_fn(cls, batch):
        keys = ['X', 'B', 'A', 'atom_positions', 'block_lengths', 'segment_ids']
        types = [torch.float, torch.long, torch.long, torch.long, torch.long, torch.long]
        res = {}
        for key, _type in zip(keys, types):
            val = []
            for item in batch:
                val.append(torch.tensor(item[key], dtype=_type))
            res[key] = torch.cat(val, dim=0)
        lengths = [len(item['B']) for item in batch]
        res['lengths'] = torch.tensor(lengths, dtype=torch.long)
        return res


class ProtInterfaceDataset(PDBDataset):
    def __init__(self, data_file: str, start_idx: int = 0, num_lines: int | None = None):
        super().__init__(data_file, start_idx=start_idx, num_lines=num_lines)

        for item in self.data:
            item['prot_data'] = BlockGeoAffDataset.filter_for_segment(item['data'], 0)

    
    def __getitem__(self, idx):
        '''
        an example of the returned data
        {
            'X': [Natom, 3],
            'B': [Nblock],
            'A': [Natom],
            'atom_positions': [Natom],
            'block_lengths': [Natom]
            'segment_ids': [Nblock]
        }        
        '''
        cmplx_data = self.data[idx]['data']
        prot_data = self.data[idx]['prot_data']

        data = {
            'cmplx': cmplx_data,
            'prot': prot_data,
            'len': len(prot_data['B']) + len(cmplx_data['B']),
        }
        return data

    @classmethod
    def collate_fn(cls, batch):
        batch_prot = super().collate_fn([item['prot'] for item in batch])
        batch_cmplx = super().collate_fn([item['cmplx'] for item in batch])

        batch = {
            'prot': batch_prot,
            'cmplx': batch_cmplx,
        }
        return batch


class LabelledPDBDataset(torch.utils.data.Dataset):

    def __init__(self, data_file):
        super().__init__()
        self.data = open_data_file(data_file)
        self.indexes = [ item['id'] for item in self.data ]  # to satify the requirements of inference.py

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        '''
        an example of the returned data
        {
            'X': [Natom, 3],
            'B': [Nblock],
            'A': [Natom],
            'atom_positions': [Natom],
            'block_lengths': [Natom]
            'segment_ids': [Nblock]
            'label': [1]
        }        
        '''
        item = self.data[idx]
        data = item['data']
        if "label" in item.keys():
            data["label"] = item["label"]
        else:
            data['label'] = item['affinity']['neglog_aff']
        return data

    @classmethod
    def collate_fn(cls, batch):
        keys = ['X', 'B', 'A', 'atom_positions', 'block_lengths', 'segment_ids', 'label']
        types = [torch.float, torch.long, torch.long, torch.long, torch.long, torch.long, torch.float]
        res = {}
        for key, _type in zip(keys, types):
            val = []
            for item in batch:
                val.append(torch.tensor(item[key], dtype=_type))
            if key == 'label':
                res[key] = torch.tensor(val, dtype=_type)
            else:
                res[key] = torch.cat(val, dim=0)
        lengths = [len(item['B']) for item in batch]
        res['lengths'] = torch.tensor(lengths, dtype=torch.long)
        return res


class MultiClassLabelledPDBDataset(torch.utils.data.Dataset):

    def __init__(self, data_file):
        super().__init__()
        self.data = open_data_file(data_file)
        self.indexes = [ item['id'] for item in self.data ]  # to satify the requirements of inference.py

    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        '''
        an example of the returned data
        {
            'X': [Natom, 3],
            'B': [Nblock],
            'A': [Natom],
            'atom_positions': [Natom],
            'block_lengths': [Natom]
            'segment_ids': [Nblock]
            'label': [1]
        }        
        '''
        item = self.data[idx]
        data = item['data']
        data["label"] = item["label"]

        return data

    @classmethod
    def collate_fn(cls, batch):
        keys = ['X', 'B', 'A', 'atom_positions', 'block_lengths', 'segment_ids', 'label']
        types = [torch.float, torch.long, torch.long, torch.long, torch.long, torch.long, torch.long]
        res = {}
        for key, _type in zip(keys, types):
            val = []
            for item in batch:
                val.append(torch.tensor(item[key], dtype=_type))
            if key == 'label':
                res[key] = torch.tensor(val, dtype=_type)
            else:
                res[key] = torch.cat(val, dim=0)
        lengths = [len(item['B']) for item in batch]
        res['lengths'] = torch.tensor(lengths, dtype=torch.long)
        return res


class MixDatasetWrapper(torch.utils.data.Dataset):
    def __init__(self, *datasets) -> None:
        super().__init__()
        self.datasets = datasets
        self.cum_len = []
        self.total_len = 0
        for dataset in datasets:
            self.total_len += len(dataset)
            self.cum_len.append(self.total_len)
        self.collate_fn = self.datasets[0].collate_fn
    
    def __len__(self):
        return self.total_len
    
    def __getitem__(self, idx):
        last_cum_len = 0
        for i, cum_len in enumerate(self.cum_len):
            if idx < cum_len:
                return self.datasets[i].__getitem__(idx - last_cum_len)
            last_cum_len = cum_len
        return None

    def _get_raw_item(self, idx):
        last_cum_len = 0
        for i, cum_len in enumerate(self.cum_len):
            if idx < cum_len:
                return self.datasets[i]._get_raw_item(idx - last_cum_len)
            last_cum_len = cum_len
        return None



class DynamicBatchWrapper(torch.utils.data.Dataset):
    def __init__(self, dataset, max_n_vertex_per_batch, max_n_vertex_per_item=None, shuffle=True) -> None:
        super().__init__()
        self.dataset = dataset
        self.indexes = [i for i in range(len(dataset))]
        self.max_n_vertex_per_batch = max_n_vertex_per_batch
        if max_n_vertex_per_item is None:
            max_n_vertex_per_item = max_n_vertex_per_batch
        self.max_n_vertex_per_item = max_n_vertex_per_item
        self.total_size = None
        self.batch_indexes = []
        self.shuffle = shuffle
        self._form_batch()

    ########## overload with your criterion ##########
    def _form_batch(self):
        if self.shuffle:
            np.random.shuffle(self.indexes)
        last_batch_indexes = self.batch_indexes
        self.batch_indexes = []

        cur_vertex_cnt = 0
        batch = []
        
        num_too_large = 0

        for i in self.indexes:
            data = self.dataset[i]
            item_len = len(data['B']) if 'B' in data else data['len']
            if item_len > self.max_n_vertex_per_item:
                num_too_large += 1
                continue
            cur_vertex_cnt += item_len
            if cur_vertex_cnt > self.max_n_vertex_per_batch:
                self.batch_indexes.append(batch)
                batch = []
                cur_vertex_cnt = item_len
            batch.append(i)
        self.batch_indexes.append(batch)
        print_log(f'Number of items too large: {num_too_large} out of {len(self.indexes)}. Remaining: {len(self.indexes) - num_too_large} batches. Created {len(self.batch_indexes)} batches.')

        if self.total_size is None:
            self.total_size = len(self.batch_indexes)
        else:
            # control the lengths of the dataset, otherwise the dataloader will raise error
            if len(self.batch_indexes) < self.total_size:
                num_add = self.total_size - len(self.batch_indexes)
                self.batch_indexes = self.batch_indexes + last_batch_indexes[:num_add]
            else:
                self.batch_indexes = self.batch_indexes[:self.total_size]

    def __len__(self):
        return len(self.batch_indexes)
    
    def __getitem__(self, idx):
        return [self.dataset[i] for i in self.batch_indexes[idx]]
    
    def collate_fn(self, batched_batch):
        batch = []
        for minibatch in batched_batch:
            batch.extend(minibatch)
        return self.dataset.collate_fn(batch)

class BalancedDynamicBatchWrapper(DynamicBatchWrapper):
    def __init__(self, dataset, max_n_vertex_per_batch, max_n_vertex_per_item=None, shuffle=True) -> None:
        self.labels = Counter([item['label'] for item in dataset])
        self.sampling_weights = [1 / self.labels[item['label']] for item in dataset]
        self.sampling_weights = np.array(self.sampling_weights) / sum(self.sampling_weights)
        super().__init__(dataset, max_n_vertex_per_batch, max_n_vertex_per_item, shuffle)
        self._form_batch()
    
    def _form_batch(self):
        chosen_indexes = np.random.choice(len(self.dataset), size=len(self.dataset), replace=True, p=self.sampling_weights)
        self.indexes = chosen_indexes
        super()._form_batch()
        print(f'Number of items in each class: {Counter([self.dataset[i]["label"] for i in self.indexes])}')

class PretrainBalancedDynamicBatchWrapper(DynamicBatchWrapper):
    def __init__(self, dataset, max_n_vertex_per_batch, max_n_vertex_per_item=None, shuffle=True) -> None:
        self.dataset_labels = [dataset._get_raw_item(idx)['modality'] for idx in range(len(dataset))]
        self.labels = Counter(self.dataset_labels)
        self.sampling_weights = [1 / np.log(self.labels[label]) for label in self.dataset_labels] # downweight the large classes by log
        self.sampling_weights = np.array(self.sampling_weights) / sum(self.sampling_weights)
        super().__init__(dataset, max_n_vertex_per_batch, max_n_vertex_per_item, shuffle)
        self._form_batch()
    
    def _form_batch(self):
        chosen_indexes = np.random.choice(len(self.dataset), size=len(self.dataset), replace=True, p=self.sampling_weights)
        self.indexes = chosen_indexes
        super()._form_batch()
        print(f'Number of items in each class: {Counter([self.dataset_labels[i] for i in self.indexes])}. Original class distribution: {self.labels}')


def parse():
    parser = argparse.ArgumentParser(description='Process data')
    parser.add_argument('--dataset', type=str, required=True, help='dataset')
    parser.add_argument('--database', type=str, default=None, help='directory of pdb data')
    return parser.parse_args()

def dataset_to_compressed_jsonl(dataset, output_file):
    with gzip.open(output_file, 'wb', compresslevel=6) as f:
        for item in dataset:
            if 'block_to_pdb_indexes' in item:
                item['block_to_pdb_indexes'] = {str(k): v for k, v in item['block_to_pdb_indexes'].items()}
            f.write(orjson.dumps(item) + b"\n")

def compressed_jsonl_to_dataset(input_file, start_idx: int = 0, num_lines: int | None = None):
    """Load JSONL.gz file with optional slicing for memory efficiency.
    
    Args:
        input_file: Path to the compressed JSONL file
        start_idx: Starting line index (0-based)
        num_lines: Maximum number of lines to load (None = all remaining lines)
    """
    dataset = []
    current_idx = 0
    end_idx = start_idx + num_lines if num_lines is not None else None
    
    with gzip.open(input_file, 'rb') as f:
        for line in f:
            if current_idx < start_idx:
                current_idx += 1
                continue
            if end_idx is not None and current_idx >= end_idx:
                break
            
            item = orjson.loads(line)
            if 'block_to_pdb_indexes' in item:
                item['block_to_pdb_indexes'] = {int(k): v for k, v in item['block_to_pdb_indexes'].items()}
            dataset.append(item)
            current_idx += 1
    
    return dataset

def open_data_file(data_file, start_idx: int = 0, num_lines: int | None = None):
    if data_file.endswith('.jsonl.gz'):
        return compressed_jsonl_to_dataset(data_file, start_idx=start_idx, num_lines=num_lines)
    elif data_file.endswith('.pkl'):
        with open(data_file, 'rb') as f:
            return pickle.load(f)
    else:
        raise ValueError('Unknown file format')
 

if __name__ == '__main__':
    args = parse()
    dataset = BlockGeoAffDataset(args.dataset)
    print(len(dataset))
    length = [len(item['B']) for item in dataset]
    print(f'interface length: min {min(length)}, max {max(length)}, mean {sum(length) / len(length)}')
    atom_length = [len(item['A']) for item in dataset]
    print(f'atom number: min {min(atom_length)}, max {max(atom_length)}, mean {sum(atom_length) / len(atom_length)}')

    item = dataset[0]
    print(item)