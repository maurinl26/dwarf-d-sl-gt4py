from typing import Tuple
from config import Config
from sl_dace.dims import I, J, K

class SLDriver:

    def __init__(self, config: Config, grid: Tuple[int]):
        self.config = config
        self.grid = grid

        # sl init

        # sl xy

    @property
    def _input_fields(self):
        return {
            "vx": {"grid": (I,J,K)}
        }

    # todo : implement call
    def __call__(self, state: dict):
        ...

