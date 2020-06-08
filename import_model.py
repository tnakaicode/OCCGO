import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import time
import os

from src.model import Model

obj = Model(cfgfile="./cfg/model.json")
print(obj.cfg["name"])
