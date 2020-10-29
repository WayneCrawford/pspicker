"""
Routines that process/interpret traces or lists of traces
"""
from .select_traces import select_traces
from .smooth_filter import smooth_filter

__all__ = ['select_traces', 'smooth_filter']
