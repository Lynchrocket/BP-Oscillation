import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from sympy.abc import x, y


class CSS:
    """CSS quantum code wrapper used by the decoder stack.

    The class keeps the representation intentionally lightweight: only the
    stabiliser matrices and logical operators are stored, and all inputs are
    normalised to binary `np.uint8` arrays.
    """

    def __init__(self, H_Z, H_X, L_Z, L_X):
        self.H_X = np.asarray(H_X, dtype=np.uint8) % 2
        self.H_Z = np.asarray(H_Z, dtype=np.uint8) % 2
        self.L_X = np.asarray(L_X, dtype=np.uint8) % 2
        self.L_Z = np.asarray(L_Z, dtype=np.uint8) % 2
        self.k, self.n = self.L_X.shape
        self.c_x = self.H_X.shape[0]
        self.c_z = self.H_Z.shape[0]

    @classmethod
    def from_qldpc(cls, code):
        """Create a CSS wrapper from a qldpc code object."""
        import qldpc
        H_Z = np.asarray(code.matrix_z, dtype=np.uint8)
        H_X = np.asarray(code.matrix_x, dtype=np.uint8)
        L_X = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.X), dtype=np.uint8)
        L_Z = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.Z), dtype=np.uint8)
        return cls(H_Z, H_X, L_Z, L_X)

    def __repr__(self):
        return f"CSS code: [[{self.n}, {self.k}, d]]"

    def draw(self):
        CSS.draw_css_tanner(self.H_X, self.H_Z, self.L_X, self.L_Z)

    def tanner(self):
        self.graph, _, _ = CSS.tanner_graph(self.H_X, self.H_Z, self.L_X, self.L_Z)
        return self.graph

    def get_matrices(self):
        return self.H_X, self.H_Z, self.L_X, self.L_Z

    @staticmethod
    def tanner_graph(
        H_X: np.ndarray,
        H_Z: np.ndarray,
        logical_Xs: np.ndarray | None = None,
        logical_Zs: np.ndarray | None = None,
    ):
        """Build a Tanner graph containing checks, data nodes and logicals."""
        n = H_X.shape[1]
        m_x, m_z = H_X.shape[0], H_Z.shape[0]
        l_x = 0 if logical_Xs is None else logical_Xs.shape[0]
        l_z = 0 if logical_Zs is None else logical_Zs.shape[0]
        G = nx.Graph()

        edge = {
            "cz": ("salmon", "dotted"),
            "cx": ("lightgreen", "dashed"),
            "lz": ("red", "solid"),
            "lx": ("green", "solid"),
        }

        pos = {}
        edge_styles = {}

        for i in range(n):
            G.add_node(f"d{i}", type="data")
            pos[f"d{i}"] = (1, n / 2 - i)

        for i in range(m_z):
            G.add_node(f"z{i}", type="z_check")
            pos[f"z{i}"] = (0, (m_z + l_z) / 2 - i)
            for j in range(n):
                if H_Z[i, j]:
                    G.add_edge(f"z{i}", f"d{j}")
                    edge_styles[(f"z{i}", f"d{j}")] = edge["cz"]

        for i in range(m_x):
            G.add_node(f"x{i}", type="x_check")
            pos[f"x{i}"] = (2, (m_x + l_x) / 2 - i)
            for j in range(n):
                if H_X[i, j]:
                    G.add_edge(f"x{i}", f"d{j}")
                    edge_styles[(f"x{i}", f"d{j}")] = edge["cx"]

        if logical_Zs is not None:
            for i, logical in enumerate(logical_Zs):
                node_name = f"Lz{i}"
                G.add_node(node_name, type="logical_z")
                pos[node_name] = (0, (m_z + l_z) / 2 - (m_z + 1 + i))
                for j in range(n):
                    if logical[j]:
                        G.add_edge(node_name, f"d{j}")
                        edge_styles[(node_name, f"d{j}")] = edge["lz"]

        if logical_Xs is not None:
            for i, logical in enumerate(logical_Xs):
                node_name = f"Lx{i}"
                G.add_node(node_name, type="logical_x")
                pos[node_name] = (2, (m_x + l_x) / 2 - (m_x + 1 + i))
                for j in range(n):
                    if logical[j]:
                        G.add_edge(node_name, f"d{j}")
                        edge_styles[(node_name, f"d{j}")] = edge["lx"]

        return G, edge_styles, pos

    @staticmethod
    def draw_css_tanner(
        H_X: np.ndarray,
        H_Z: np.ndarray,
        logical_Xs: np.ndarray | None = None,
        logical_Zs: np.ndarray | None = None,
    ):
        G, edge_styles, pos = CSS.tanner_graph(H_X, H_Z, logical_Xs, logical_Zs)
        n = H_X.shape[1]

        node_colors = {
            "data": "skyblue",
            "z_check": "salmon",
            "x_check": "lightgreen",
            "logical_z": "red",
            "logical_x": "green",
        }
        node_shapes = {
            "data": "o",
            "z_check": "s",
            "x_check": "s",
            "logical_z": "o",
            "logical_x": "o",
        }

        fig, ax = plt.subplots(figsize=(10, max(6, n / 2)))

        for node_type, shape in node_shapes.items():
            nodes = [node for node in G.nodes if G.nodes[node].get("type") == node_type]
            nx.draw_networkx_nodes(
                G,
                pos,
                nodelist=nodes,
                node_shape=shape,
                node_color=node_colors[node_type],
                label=node_type,
                ax=ax,
            )

        for (u, v), (color, style) in edge_styles.items():
            nx.draw_networkx_edges(
                G,
                pos,
                edgelist=[(u, v)],
                edge_color=color,
                style=style,
                ax=ax,
            )

        nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
        ax.set_axis_off()
        plt.tight_layout()
        plt.show()


class IcebergCode(CSS):
    """[[n, n-2, 2]] iceberg code, valid only for even n"""

    def __init__(self, n):
        assert n % 2 == 0
        H_Z = np.asarray([[1 for _ in range(n)]], dtype=np.uint8)
        H_X = np.asarray([[1 for _ in range(n)]], dtype=np.uint8)
        L_X = np.asarray(
            [[1] + [1 if j == i else 0 for j in range(n - 1)] for i in range(n - 2)],
            dtype=np.uint8,
        )
        L_Z = np.asarray(
            [[1 if j == i else 0 for j in range(n - 1)] + [1] for i in range(1, n - 1)],
            dtype=np.uint8,
        )
        super().__init__(H_Z, H_X, L_Z, L_X)


class SteaneCode(CSS):
    """[[7, 1, 3]] Steane code"""

    def __init__(self):
        H = [
            [0, 1, 1, 1, 1, 0, 0],
            [1, 0, 1, 1, 0, 1, 0],
            [1, 1, 0, 1, 0, 0, 1],
        ]
        H_X = np.asarray(H, dtype=np.uint8)
        H_Z = np.asarray(H, dtype=np.uint8)
        L_Z = np.asarray([[1 for _ in range(7)]], dtype=np.uint8)
        L_X = np.asarray([[1 for _ in range(7)]], dtype=np.uint8)
        super().__init__(H_Z, H_X, L_Z, L_X)

class SurfaceCode(CSS):
    def __init__(self, d):
        import qldpc
        code = qldpc.codes.SurfaceCode(d)
        H_Z = np.asarray(code.matrix_z)
        H_X = np.asarray(code.matrix_x)
        L_X = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.X))
        L_Z = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.Z))
        super(SurfaceCode, self).__init__(H_Z, H_X, L_Z, L_X)
        

class ToricCode(CSS):
    def __init__(self, d):
        import qldpc
        rep_code = qldpc.codes.RingCode(d)
        code = qldpc.codes.HGPCode(rep_code, rep_code)
        H_Z = np.asarray(code.matrix_z)
        H_X = np.asarray(code.matrix_x)
        L_X = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.X))
        L_Z = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.Z))
        super(ToricCode, self).__init__(H_Z, H_X, L_Z, L_X)


class BivariateBicycleCode(CSS):
    def __init__(self, orders, poly_a, poly_b):
        import qldpc
        code = qldpc.codes.BBCode(orders, poly_a, poly_b)
        H_Z = np.asarray(code.matrix_z)
        H_X = np.asarray(code.matrix_x)
        L_X = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.X))
        L_Z = np.asarray(code.get_logical_ops(qldpc.objects.Pauli.Z))
        super(BivariateBicycleCode, self).__init__(H_Z, H_X, L_Z, L_X)

class GrossCode(BivariateBicycleCode):
    def __init__(self):
        orders = {x: 12, y: 6}
        poly_a = x**3 + y + y**2
        poly_b = y**3 + x + x**2
        super(GrossCode, self).__init__(orders, poly_a, poly_b)
