from dataclasses import dataclass

@dataclass
class Point:
    x: float = 0
    y: float = 0
    def __repr__(self) -> str:
        return "(%f, %f)" % (self.x, self.y)