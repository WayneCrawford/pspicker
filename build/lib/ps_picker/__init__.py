"""
Routines to pick P and S waves using kurtosis extrema and to
"associate" them using clustering

Can also use the Kurtosis class independently
"""
from .ps_picker import PSPicker
from .kurtosis import Kurtosis

__all__ = ['PSPicker', 'Kurtosis']
