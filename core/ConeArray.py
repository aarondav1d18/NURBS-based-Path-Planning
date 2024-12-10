from dataclasses import dataclass
from typing import List, Tuple
from core.Cone import Cone


@dataclass
class ConeArrayAndGraph:
	cone_stack: List[Cone]
	cone_graph: List[bool]