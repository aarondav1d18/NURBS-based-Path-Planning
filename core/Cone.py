from dataclasses import dataclass
from enum import Enum
from typing import List
from core.Point import Point


class ConeColour(Enum) :
	YELLOW = 1,
	BLUE = 2,
	UNKNOWN = 3
	ORANGE = 4
	
@dataclass
class Cone:
	pos: Point = None
	colour: ConeColour = ConeColour.UNKNOWN
	covariance: List[List[float]] = None

	def __repr__(self) -> str:
		# return "(%f, %f, %s)" % (self.x, self.y, self.colour)
		return "(%f, %f)" % (self.pos.x, self.pos.y)
	
	def __hash__(self):
		# Hash based on unique position of the cone
		return hash((self.pos.x, self.pos.y))

	def __eq__(self, other):
		# Check equality based on position of the cone
		return isinstance(other, Cone) and self.pos.x == other.pos.x and self.pos.y == other.pos.y