import sys
import types
import functools

from typing import Sequence, Optional, Any


def flatten(seq: Sequence[Optional[Any]]):
    return [item for subseq in seq for item in subseq]


def is_callable(func):
    """Check if the object is callable."""
    FUNCTYPES = (types.FunctionType, types.MethodType, functools.partial)
    return func and (isinstance(func, FUNCTYPES) or callable(func))
