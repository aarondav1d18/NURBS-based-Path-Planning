from dataclasses import dataclass

@dataclass
class PointVelocity:
    x: float = 0
    y: float = 0
    vel: float = 0
    def __repr__(self) -> str:
        return "(%f, %f)" % (self.x, self.y)