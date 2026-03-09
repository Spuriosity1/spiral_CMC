import h5py
import numpy as np

from fury import window, actor
from matplotlib import cm
from matplotlib import colormaps
from matplotlib.colors import Normalize


def add_reference_axes(scene, length=1.0, linewidth=3):
    """
    Add XYZ reference axes to a FURY scene.

    X = Red
    Y = Green
    Z = Blue
    """

    origin = np.array([[0, 0, 0]])

    # Unit directions
    directions = np.array([
        [length, 0, 0],   # X
        [0, length, 0],   # Y
        [0, 0, length],   # Z
    ])

    # Colors (RGB)
    colors = np.array([
        [1, 0, 0],  # X = red
        [0, 1, 0],  # Y = green
        [0, 0, 1],  # Z = blue
    ])

    starts = np.repeat(origin, 3, axis=0)

    axes = actor.arrow(
        starts,
        directions,
        colors=colors,
        scales=linewidth
    )

    scene.add(axes)

class SpinHDF5Parser:
    """
    Parser and visualizer for spin HDF5 files using FURY.
    """

    def __init__(self, filename: str):
        self.filename = filename

        self.positions = None
        self.orientations = None

    # ------------------------------------------------------
    # Public API
    # ------------------------------------------------------

    def load(self):
        """Load and validate HDF5 data."""
        self._read_file()
        self._validate()
        return self

    def visualise(
        self,
        scale: float = 1.0,
        normalize: bool = True,
        cmap: str = "viridis",
        bg_color=(0, 0, 0),
        window_size=(1000, 800),
    ):
        """
        Visualize spins in 3D.

        Parameters
        ----------
        scale : float
            Vector length scaling
        normalize : bool
            Normalize spin vectors
        cmap : str
            Matplotlib colormap
        """

        pos = self.positions.copy()
        vec = self.orientations.copy()

        if normalize:
            vec = self._normalize(vec)

        vec *= scale

        # Compute scalar field (magnitudes)
        mags = np.linalg.norm(vec, axis=1)

        # Map scalars -> RGB
        colors = self._apply_colormap(mags, cmap)

        scene = window.Scene()
        scene.background(bg_color)

        arrows = actor.arrow(
            pos,
            vec,
            colors=colors,
            scales=1
        )

        scene.add(arrows)

        add_reference_axes(scene)

        window.show(
            scene,
            size=window_size,
            reset_camera=True
        )

    # ------------------------------------------------------
    # Internal
    # ------------------------------------------------------

    def _read_file(self):

        try:
            with h5py.File(self.filename, "r") as f:

                if "spin_pos" not in f:
                    raise KeyError("Missing spin_pos dataset")

                if "spin_orientation" not in f:
                    raise KeyError("Missing spin_orientation dataset")

                self.positions = np.asarray(f["spin_pos"])
                self.orientations = np.asarray(f["spin_orientation"])

        except OSError as e:
            raise RuntimeError(f"HDF5 open failed: {e}")

    def _validate(self):

        p = self.positions
        o = self.orientations

        if p.ndim != 2 or o.ndim != 2:
            raise ValueError("Datasets must be 2D")

        if p.shape[1] != 3:
            raise ValueError("spin_pos must be (N,3)")

        if o.shape[1] != 3:
            raise ValueError("spin_orientation must be (N,3)")

        if p.shape[0] != o.shape[0]:
            raise ValueError("Size mismatch")

        if not np.isfinite(p).all():
            raise ValueError("Invalid spin_pos values")

        if not np.isfinite(o).all():
            raise ValueError("Invalid spin_orientation values")

    def _normalize(self, vec):

        norms = np.linalg.norm(vec, axis=1, keepdims=True)
        norms[norms == 0] = 1.0

        return vec / norms

    def _apply_colormap(self, scalars, cmap_name):
        """
        Convert scalar values to RGB colors.
        """

        norm = Normalize(
            vmin=np.min(scalars),
            vmax=np.max(scalars)
        )

#        cmap = colormaps[cmap_name]

#        rgba = cmap(norm(scalars))


        # Drop alpha channel → (N,3)
        rgb = np.array([1.,1.,1.])

        return rgb.astype(np.float32)





# ------------------------------------------------------
# CLI
# ------------------------------------------------------

def main():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("file")
    parser.add_argument("--scale", type=float, default=1.0)
    parser.add_argument("--no-normalize", action="store_true")
    parser.add_argument("--cmap", default="viridis")

    args = parser.parse_args()

    vis = SpinHDF5Parser(args.file).load()

    vis.visualise(
        scale=args.scale,
        normalize=not args.no_normalize,
        cmap=args.cmap
    )


if __name__ == "__main__":
    main()

