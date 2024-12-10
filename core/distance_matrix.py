import math
from core.Point import Point
from core.Cone import Cone

def euclidean_distance_between_cones(cone_a: Cone, cone_b: Cone) -> float :
	# just adding checks for if either are cones it will convert them to points
	# if type(cone_a) is Cone:
	# 	dummy = cone_a
	# 	cone_a = Point()
	# 	cone_a.x = dummy.pos.x
	# 	cone_a.y = dummy.pos.y
	# elif type(cone_b) is Cone:
	# 	dummy = cone_b
	# 	cone_b = Point()
	# 	cone_b.x = dummy.pos.x
	# 	cone_b.y = dummy.pos.y
	if type(cone_a) is Point:
		dummy = cone_a
		cone_a = Cone()
		cone_a.pos = dummy
	elif type(cone_b) is Point:
		dummy = cone_b
		cone_b = Cone()
		cone_b.pos = dummy
	return math.hypot(cone_a.pos.x - cone_b.pos.x, cone_a.pos.y - cone_b.pos.y) 
