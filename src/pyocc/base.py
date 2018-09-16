from abc import abstractmethod
from collections import MutableMapping
import os

class Base(MutableMapping, object):
    """
    Notes
    -----
    When properly defined in inherited functions, this class should behave like
    a dictionary.

    As this class inherits from MutableMapping, any class inherting from
    AirconicsBase must also define the abstract methods of Mutable mapping,
    i.e. __setitem__, __getitem__, __len__, __iter__, __delitem__
    """
    @abstractmethod
    def __init__(self, *args, **kwargs):
        pass

    @abstractmethod
    def __str__(self, *args, **kwargs):
        pass

    @abstractmethod
    def Display(self, *args, **kwargs):
        pass

    @abstractmethod
    def Write(self, *args, **kwargs):
        pass
    
    @abstractmethod
    def Build(self, *args, **kwargs):
        pass