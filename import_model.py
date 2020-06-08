import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import time
import os

from src.model import model, model_base

obj = model(cfgfile="./cfg/model.json")
print(obj.cfg["name"])

md0 = obj.set_model(name="surf0")
md1 = obj.set_model(name="surf1")
md2 = obj.set_model(name="surf")

print(md0.axs.Location())
print(md2.axs.Location())
