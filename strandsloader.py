import torch
import numpy as np
import random 

class StrandLoader():
    def __init__(self, data):
        self.data = data
        self.data_len = len(data)

    def __iter__(self):
        raise NotImplementedError("SubClass should implement this method .")

class SequentialStrandLoader(StrandLoader):

    def __init__(self, data):
        super().__init__(data)
        self.index = 0

    def __iter__(self):
        return self
    
    def __next__(self):
        if self.index < self.data_len:
            item = self.data[self.index]
            self.index += 1
            return item
        else :
            raise StopIteration

class ReverseStrandLoader(StrandLoader):
    def __init__(self, data):
        super().__init__(data)
        self.index = self.data_len - 1
    
    def __iter__(self):
        return self
        
    def __next__(self):
        if self.index >= 0:
            item = self.data[self.index]
            self.index -= 1
            return item
        else:
            raise StopIteration
           
class PermutationsStrandLoader(StrandLoader):
    def __init__(self, data, seed):
        super().__init__(data)
        self.seed = seed
        

    def __iter__(self):
        random.seed(self.seed)
        self.permutations = np.random.permutation(len(self.data))
        self.iter_index = 0
        return self
    
    def __next__(self):
        if self.iter_index < len(self.permutations):
            item = self.data[self.permutations[self.iter_index]]
            self.iter_index += 1
            return item

class RandomStrandLoader(StrandLoader):
    def __init__(self, data, seed):
        super().__init__(data)
        self.indices = list(range(len(data)))
        self.seed = seed
    
    def __iter__(self):
        random.seed(self.seed)
        random.shuffle(self.indices)  # Shuffle indices for new iteration
        self.iter_index = 0
        return self

    def __next__(self):
        if self.iter_index < len(self.indices):
            item = self.data[self.indices[self.iter_index]]
            self.iter_index += 1
            return item
        else:
            raise StopIteration