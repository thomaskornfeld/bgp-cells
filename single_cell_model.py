import numpy as np

class Environment:
    """
    Represents the discrete 2D concentration field c(x,y).
    """
    def __init__(self, concentration_matrix, dx=1.0):
        self.c_field = np.array(concentration_matrix)
        self.nx, self.ny = self.c_field.shape
        self.dx = dx # Spatial scale of a single discrete grid cell

    def get_discrete_pos(self, x, y):
        """Maps continuous (x,y) to discrete grid indices."""
        # Round continuous position to nearest discrete grid index
        i = int(np.clip(np.round(x / self.dx), 0, self.nx - 1))
        j = int(np.clip(np.round(y / self.dx), 0, self.ny - 1))
        return i, j

    def get_concentration(self, x, y):
        """Returns the concentration at the cell's current position."""
        i, j = self.get_discrete_pos(x, y)
        return self.c_field[i, j]


class Cell:
    """
    Simulates a cell with spatial Brownian motion and 
    stochastic gene expression driven by an OU process.
    """
    def __init__(self, x0, y0, g0, q0, beta, gamma, D_q, D_x, mu_func):
        # State variables
        self.x = x0
        self.y = y0
        self.g = g0
        self.q = q0
        
        # Parameters
        self.beta = beta       # Gene degradation rate
        self.gamma = gamma     # OU process relaxation rate
        self.D_q = D_q         # Diffusion coefficient for q noise
        self.D_x = D_x         # Spatial diffusion coefficient
        self.mu_func = mu_func # Function mapping concentration c to mean target mu(c)
        
        # Arrays to store history for analysis
        self.history = {
            'x': [self.x],
            'y': [self.y],
            'g': [self.g],
            'q': [self.q]
        }

    def step(self, dt, env):
        """Advances the cell's state by one time step dt."""
        # 1. Spatial Brownian Motion
        # dX = sqrt(2 * D_x) * dW
        self.x += np.sqrt(2 * self.D_x * dt) * np.random.randn()
        self.y += np.sqrt(2 * self.D_x * dt) * np.random.randn()

        # 2. Get local concentration to determine target mean mu
        c_local = env.get_concentration(self.x, self.y)
        mu_target = self.mu_func(c_local)

        # 3. Update production rate q (Ornstein-Uhlenbeck process via Euler-Maruyama)
        # dq = gamma * (mu(c) - q) * dt + sqrt(2 * D_q) * dW
        dq = self.gamma * (mu_target - self.q) * dt + np.sqrt(2 * self.D_q * dt) * np.random.randn()
        self.q += dq

        # 4. Update gene product g (Deterministic Euler step based on q)
        # dg = (q - beta * g) * dt
        dg = (self.q - self.beta * self.g) * dt
        self.g += dg

        # Record history
        self.history['x'].append(self.x)
        self.history['y'].append(self.y)
        self.history['g'].append(self.g)
        self.history['q'].append(self.q)