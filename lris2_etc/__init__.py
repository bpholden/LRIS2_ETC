"""Exposure time calculator for LRIS-2"""
import os

# templates, filters, sky and throughput files ship inside the package
DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
