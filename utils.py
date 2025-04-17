from numpy import pi


def unify_angle(angle):
    """Normalize angle to [-180, 180]"""
    while angle > 180:
        angle -= 360
    while angle < -180:
        angle += 360
    return angle


def calculate_dihedral_angle(self, atom1, atom2, atom3, atom4):
    self.find_res_in_structure()
    angle = calculate_dihedral_angle(
        atom1.get_vector(),
        atom2.get_vector(),
        atom3.get_vector(),
        atom4.get_vector()
    ) * 180 / pi
    return angle



