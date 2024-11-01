import os.path as osp
import torch
from torch_geometric.data import InMemoryDataset

class MyDataset(InMemoryDataset):
    """
    A PyTorch InMemoryDataset to build multi-view dataset through graph data augmentation
    """
    def __init__(self, root="../datasets/", file=None
                 , transform=None, pre_transform=None, pre_filter=None):
        super().__init__(root=root+file, transform=transform, pre_transform=pre_transform, pre_filter=pre_filter)
        self.data, self.slices = torch.load(self.processed_paths[0])

    def process(self):
        print("Processing full batch data")
        data = torch.load(osp.join(self.raw_dir, self.raw_file_names[0]))
        data_list = [data]
        if self.pre_filter is not None:
            data_list = [data for data in data_list if self.pre_filter(data)]
        if self.pre_transform is not None:
            data_list = [self.pre_transform(data) for data in data_list]
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])

    def download(self):
        pass

    @property
    def raw_dir(self) -> str:
        return osp.join(self.root, 'raw')
    
    @property
    def raw_file_names(self):
        return ['raw.pt']
    
    @property
    def processed_dir(self) -> str:
        return osp.join(self.root, 'processed')
    
    @property
    def processed_file_names(self):
        return ['processed.pt']