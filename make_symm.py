import numpy as np
from solid import scad_render_to_file, sphere, translate


class PointGroup:
    def rot_mat(theta, axis):
        axis = np.mat([axis])
        return (np.cos(theta) * np.identity(3)
            + (1 - np.cos(theta)) * axis.T * axis
            + np.sin(theta) * np.cross(axis, np.identity(3)))
    # 主軸はzとする
    # C2はxとする
    def E():
        return [np.identity(3)]
    def C(axis, n):
        theta = 2 * np.pi / n
        return list(map(
            lambda i: PointGroup.rot_mat(theta * i, axis),
            range(n)
        ))
    def Ci():
        return -np.identity(3)
    def sigma(axis):
        axis = np.mat([axis])
        return [np.identity(3) - 2 * axis.T * axis]
    def sigma_h():
        return [PointGroup.sigma(np.array([0, 0, 1]))]
    def sigma_v(n):
        y_axis = np.mat([[0, 1, 0]]).T
        sigma = []
        for c in PointGroup.C([0, 0, 1], n):
            sigma.extend(PointGroup.sigma(np.ravel(c * y_axis)))
        return sigma
    def sigma_d(n):
        theta2 = np.pi / n / 2
        rot = PointGroup.rot_mat(theta2, [0, 0, 1])
        rot_inv = rot ** -1
        return list(map(lambda m: rot * m * rot_inv, PointGroup.sigma_v(n)))
    def S(n):
        axis = np.array([0, 0, 1])
        sigma = PointGroup.sigma_h()[0]
        rot = PointGroup.C(axis, n)[1]
        S = sigma * rot
        return list(map(lambda i: S ** i, range(n)))
    def D(n):
        axis = np.array([0, 0, 1])
        sub_axis = [1, 0, 0]
        rot = PointGroup.rot_mat(np.pi / n, axis)
        D = PointGroup.C(axis, n)
        sub_axis = np.mat([sub_axis]).T
        D.extend(sum(map(
            lambda i: PointGroup.C(np.ravel(rot ** i * sub_axis), 2),
            range(n)), []))
        return D
    def Dnh(n):
        D = PointGroup.D(n)
        D.extend(PointGroup.sigma_h())
        return D
    def Dnd(n):
        D = PointGroup.D(n)
        D.extend(PointGroup.sigma_d(n))
        return D



class Points:
    def __init__(self, points, mats):
        self.points = points
        self.mats = mats
    def apply_mat10(self):
        for _ in range(10):
            is_changed = self.apply_mat()
            if not is_changed:
                return False
        return True
    def apply_mat(self):
        append_points = []
        for p in self.points:
            new_points = self.apply_mat_for_p(p)
            for new_p in new_points:
                append_points.extend(new_points)
        for p in append_points:
            if not self.is_there_close_p(p):
                self.points.append(p)
        return len(append_points) != 0
    def apply_mat_for_p(self, p):
        p = np.mat([p]).T
        return list(map(lambda m: np.ravel(m * p), self.mats))

    def is_there_close_p(self, p):
        return any(map(lambda p2: all(np.isclose(p, p2, rtol=1e-3)), self.points))

    def get_spheres(self):
        return list(map(Points.put_sphere, self.points))
    def put_sphere(p):
        c = sphere(2, segments=10)
        c = translate(p)(c)
        return c

points = [[0, 0, 0], [0, 0, 10], [10, 10, 20]]
mats = PointGroup.Dnd(2)
p = Points(points, mats)
p.apply_mat10()
spheres = p.get_spheres()
mol = spheres[0]
for s in spheres:
    mol += s
scad_render_to_file(mol, 'symm_hex.scad', include_orig_code=False)
