import numpy as np

class Environment:
    """
    Represents the discrete 2D concentration field c(x,y).
    """
    def __init__(self, concentration_matrix, dx=1.0):
        self.c_field = np.array(concentration_matrix)
        self.nx, self.ny = self.c_field.shape
        self.dx = dx

    def get_discrete_pos(self, x, y):
        i = int(np.clip(np.round(x / self.dx), 0, self.nx - 1))
        j = int(np.clip(np.round(y / self.dx), 0, self.ny - 1))
        return i, j

    def get_concentration(self, x, y):
        i, j = self.get_discrete_pos(x, y)
        return self.c_field[i, j]


class Cell:
    """
    Simulates a cell with spatial Brownian motion, stochastic gene expression,
    and the ability to divide.
    """
    def __init__(self, cell_id, x0, y0, g0, q0, beta, gamma, D_q, D_x, mu_func, parent_id=None):
        self.cell_id = cell_id
        self.parent_id = parent_id
        self.x = x0
        self.y = y0
        self.g = g0
        self.q = q0
        
        self.beta = beta       
        self.gamma = gamma     
        self.D_q = D_q         
        self.D_x = D_x         
        self.mu_func = mu_func 

    def step(self, dt, env):
        # 1. Spatial Brownian Motion
        self.x += np.sqrt(2 * self.D_x * dt) * np.random.randn()
        self.y += np.sqrt(2 * self.D_x * dt) * np.random.randn()

        # 2. Local concentration and target mean
        c_local = env.get_concentration(self.x, self.y)
        mu_target = self.mu_func(c_local)

        # 3. Update production rate q (Euler-Maruyama)
        dq = self.gamma * (mu_target - self.q) * dt + np.sqrt(2 * self.D_q * dt) * np.random.randn()
        self.q += dq

        # 4. Update gene product g (Euler step)
        dg = (self.q - self.beta * self.g) * dt
        self.g += dg

    def divide(self, new_cell_id):
        """
        Creates a daughter cell inheriting the current state.
        The parent cell (self) also continues to exist.
        """
        daughter = Cell(
            cell_id=new_cell_id,
            x0=self.x, y0=self.y,
            g0=self.g, q0=self.q,
            beta=self.beta, gamma=self.gamma,
            D_q=self.D_q, D_x=self.D_x,
            mu_func=self.mu_func,
            parent_id=self.cell_id
        )
        return daughter