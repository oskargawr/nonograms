import pyswarm
import numpy as np
import sys
import os
import json
from nonograms import Nonogram
import pandas as pd

nonograms = []

def load_nonograms():
    for filename in os.listdir("nonograms"):
        if filename.endswith(".json"):
            with open(f"db/{filename}", "r") as file:
                data = json.load(file)
                nonograms.append(Nonogram(data["columns"], data["rows"], data["solution"]))


print(nonograms[4:])