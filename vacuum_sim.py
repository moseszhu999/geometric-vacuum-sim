import numpy as np
import matplotlib.pyplot as plt


class VacuumLattice:
    def __init__(self, resolution=100, scale=10.0):
        """
        Initialize the vacuum lattice (Space-Time Pixels).
        resolution: Number of pixels per dimension (N x N x N grid).
        scale: Physical size of the simulation box (in Planck lengths).
        """
        self.resolution = resolution
        self.scale = scale
        self.grid_size = scale / resolution
        
        # Create the 3D grid points
        x = np.linspace(-scale/2, scale/2, resolution)
        y = np.linspace(-scale/2, scale/2, resolution)
        z = np.linspace(-scale/2, scale/2, resolution)
        self.X, self.Y, self.Z = np.meshgrid(x, y, z, indexing='ij')
        
        # Initialize the field (Superfluid Velocity / E-field)
        self.Vx = np.zeros_like(self.X)
        self.Vy = np.zeros_like(self.Y)
        self.Vz = np.zeros_like(self.Z)
        
        print(f"Vacuum Lattice Initialized: {resolution}^3 pixels")
        print(f"Pixel Size: {self.grid_size:.4f} l_p")

    def add_vortex_ring(self, radius=2.0, circulation=1.0, core_radius=0.2):
        """
        Add a toroidal vortex ring (Electron model) to the lattice.
        radius: Major radius of the torus (R).
        circulation: Strength of the vortex (Gamma).
        core_radius: Minor radius of the torus (a).
        """
        print(f"Adding Vortex Ring: R={radius}, a={core_radius}")
        
        # Calculate distance from the z-axis (cylindrical rho)
        rho = np.sqrt(self.X**2 + self.Y**2)
        
        # Distance from the ring core center
        # The core is a circle in the xy-plane at radius R
        # d_core is the distance from any point (x,y,z) to the circle R
        d_core = np.sqrt((rho - radius)**2 + self.Z**2)
        
        # Vortex velocity field (Biot-Savart law approximation for a ring)
        # Simplified model: Tangential velocity around the core
        # v = Gamma / (2 * pi * d_core)
        # We add a core regularization to avoid infinity
        v_mag = circulation / (2 * np.pi * (d_core + 0.01))
        
        # Direction: Tangential to the minor circle (wrapping around the ring wire)
        # The flow circulates around the wire.
        # Vector from core center to point:
        # Radial vector in minor cross-section: (rho - R, z) normalized
        # Tangential vector is perpendicular to this.
        
        # Unit vectors
        # Radial direction in xy plane: (x/rho, y/rho, 0)
        ur_x = np.where(rho > 0, self.X/rho, 0)
        ur_y = np.where(rho > 0, self.Y/rho, 0)
        
        # Vector from core to point (in the poloidal plane)
        # dr = rho - R
        # dz = Z
        # We want the cross product of the toroidal direction (phi) and the radial direction?
        # No, the flow is poloidal (wrapping around the ring).
        # Wait, is the electron flow toroidal or poloidal?
        # "Tornado of Light" usually implies toroidal flow (spinning like a ring) AND poloidal (twisting).
        # Let's model a simple toroidal flow first (like a smoke ring moving).
        # Actually, a smoke ring has poloidal flow (rolling).
        # Let's assume the "flux" is the magnetic field lines, which wrap around the current loop.
        # If the electron is a current loop, the B-field wraps around it (Poloidal).
        # If the electron is a photon loop, the photon travels Toroidally.
        
        # Let's model the PHOTON FLUX traveling TOROIDALLY (around the major radius).
        # v_toroidal = c (speed of light) inside the core.
        
        # Mask for the core
        in_core = d_core < core_radius
        
        # Toroidal direction (-y, x, 0)
        # v_phi = (-y, x, 0) / rho
        v_phi_x = -ur_y
        v_phi_y = ur_x
        
        # Assign velocity only inside the core (Soliton)
        self.Vx[in_core] = v_phi_x[in_core] * circulation
        self.Vy[in_core] = v_phi_y[in_core] * circulation
        # Vz remains 0 for pure toroidal flow

    def calculate_impedance(self):
        """
        Calculate the Effective Geometric Impedance (Z_eff).
        Z_eff = Total Energy / (Total Flux)^2
        """
        # 1. Calculate Total Energy in the Lattice
        # Energy Density u ~ V^2 (Kinetic Energy of Superfluid)
        # We sum v_mag^2 over all pixels
        V_mag_sq = self.Vx**2 + self.Vy**2 + self.Vz**2
        total_energy = np.sum(V_mag_sq) * (self.grid_size**3)
        
        # 2. Calculate Total Flux (Phi) through the Ring
        # We integrate the field passing through the x-z plane (y=0) inside the ring
        # The ring is in the xy plane? No, let's check add_vortex_ring.
        # It calculates rho = sqrt(x^2 + y^2), so the ring is in the xy plane.
        # The flow wraps around the wire (poloidal).
        # So the flux passes through the hole of the donut (z-axis)?
        # No, if flow is poloidal (around the wire), the flux is "trapped" in the wire.
        # Let's calculate the flux passing through a cross-section of the wire.
        # Cross section: x-z plane at y=0, for x > R-a and x < R+a.
        
        # Slice at y index corresponding to y=0
        mid_y = self.resolution // 2
        # Flux is the flow perpendicular to the cross-section?
        # The flow is in the (x,z) plane at y=0.
        # Vx and Vz components.
        # We sum sqrt(Vx^2 + Vz^2) in the core region.
        
        flux_slice = np.sqrt(self.Vx[:, mid_y, :]**2 + self.Vz[:, mid_y, :]**2)
        # Filter for core only (approximate)
        # We can just sum the magnitude, assuming it's the flux we care about.
        total_flux = np.sum(flux_slice) * (self.grid_size**2)
        
        if total_flux == 0:
            return 0
            
        # Z_eff = E / Phi^2
        z_eff = total_energy / (total_flux**2)
        
        print(f"Total Energy: {total_energy:.4e}")
        print(f"Total Flux: {total_flux:.4e}")
        print(f"Geometric Impedance Z_eff: {z_eff:.4e}")
        
        return z_eff

    def simulate_reflection(self):
        """
        Simulate a wave hitting the vortex core to calculate Reflection Coefficient.
        Uses the refractive index profile n(r) ~ sqrt(1 + V^2).
        """
        # 1. Generate Refractive Index Map
        # n = 1 + alpha * Field_Strength (Simplified Kerr effect)
        V_mag = np.sqrt(self.Vx**2 + self.Vy**2 + self.Vz**2)
        # Normalize V_mag to have a peak of ~137 (just for testing strong field) or keep as is
        # Let's assume n follows the Born-Infeld form: n = 1 / sqrt(1 - E^2)
        # Here we treat V as E. If V < 1, n is real.
        # Let's use a simplified linear model for the simulation: n = 1 + V_mag
        n_map = 1.0 + V_mag
        
        # 2. Extract 1D profile through the core
        # Line passing through x-axis (y=0, z=0)
        mid_y = self.resolution // 2
        mid_z = self.resolution // 2
        n_profile = n_map[:, mid_y, mid_z]
        
        # 3. Calculate Reflection Coefficient (Fresnel equation for normal incidence)
        # We approximate the core as a step function or use the peak value
        n_peak = np.max(n_profile)
        n_vac = 1.0
        
        # R = ((n2 - n1) / (n2 + n1))^2
        R = ((n_peak - n_vac) / (n_peak + n_vac))**2
        
        print(f"Peak Refractive Index: {n_peak:.4f}")
        print(f"Reflection Coefficient R: {R:.4e}")
        
        return R

if __name__ == "__main__":
    sim = VacuumLattice(resolution=100, scale=10.0)
    sim.add_vortex_ring(radius=3.0, circulation=1.0, core_radius=0.5)
    
    print("-" * 30)

    
    print("-" * 30)
    z_eff = sim.calculate_impedance()
    
    print("-" * 30)
    R = sim.simulate_reflection()
    
    # Check if Z_eff ratio is close to 137 (Scaling needed)
    # This is a qualitative check for now.
