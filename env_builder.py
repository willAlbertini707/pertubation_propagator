"""
William J. Albertini
California Polytechnic Aerospace department
AERO 452-01-2238

This file is used build .env file that the rust program uses.
Makes a python interface for built rust propagator.
"""

import os, sys
import numpy as np


class EnvBuilder:

    def __init__(self, r0: np.ndarray, v0: np.ndarray, area: float, mass: float, init_jd: float, 
                cd: float = 2.2, cr: float = 1.2, mu: float = 398600.435436, body_radius: float = 6378) -> None:
        # function __init__
        # set initial system variables
        # -------------------------------
        # r0: initial position vector km
        # v0: initial velocity vector km/s
        # area: area of spacecraft (m2)
        # mass: mass (kg)
        # init_jd: initial julian date
        # cd: coefficient of drag for spacecraft
        # cr: reflectivity constant
        # mu: gravitational parameter of earth km units
        # body_radius: celestial body radius

        self._r = r0
        self._v = v0
        self._area = area
        self._cd = cd
        self._mass = mass
        self._mu = mu
        self._init_jd = init_jd
        self._cr = cr
        self._body_radius = body_radius

    def propagate(self, t0: float, tf: float, dx: float, file_path: str = "tmp.txt", add_pertubations: str = "true") -> None:
        # function propogate

        env_contents = f"""
initial_state={self._r[0]},{self._r[1]},{self._r[2]},{self._v[0]},{self._v[1]},{self._v[2]}
area={self._area}
cd={self._cd}
mass={self._mass}
mu={self._mu}
t0={t0}
tf={tf}
dx={dx}
output_file={file_path}
init_jd={self._init_jd}
cr={self._cr}
body_radius={self._body_radius}
add_pertubations={add_pertubations}
"""
        # write .env file
        self.write_env(env_contents)

        # execute rust program
        os.system("./target/release/pertubation_propagator")

        # clean .env
        # self.clean_env()


    def write_env(self, file_contents: str) -> None:
        # function write_env
        # ----------------------
        # file_contents: string to write to .env file

        if os.path.exists(".env"):
            self.clean_env()

        with open(".env", "w+") as f:
            f.write(file_contents)

        print("Environment created")

    
    def clean_env(self) -> None:
        # function clean_env
        # -------------------
        # delete .env file

        os.remove(".env")