import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt


# Constants
Ec = 32867e6 # Elastic modulus of Concrete [N/m^2 or Pa]
num_floors = 6
h = 3.2     #height of each floor
# Mass matrix (diagonal)
M = [
    [547096.5379, 0, 0, 0,0,0],
    [0, 535883.49, 0, 0,0,0 ],
    [0, 0, 525689.8101, 0,0,0 ],
    [0, 0, 0, 512692.8682,0,0 ],
    [0, 0, 0, 0,497555.2535,0 ],
    [0, 0, 0, 0,0,191512.1884 ]
    ]

# Stiffness Matrix
K = [
    [2772178314, -1147283905, 0, 0, 0,0],
    [-1147283905, 1930894248, -783610344, 0, 0,0],
    [0, -783610344, 1297737091, -514126746.7, 0,0],
    [0, 0, -514126746.7, 702271590.2, -188144843.6,0],
    [0, 0, 0, -188144843.6, 249078383.9,-60933540.34],
    [0, 0, 0, 0, -60933540.34, 60933540.34]
    ]

# Solve eigenvalue problem
eigenvalues, eigenvectors = eigh(K, M)

# Natural frequencies and periods
print(eigenvalues)
omega_n = np.sqrt(np.abs(eigenvalues))  # rad/sec
print("Natural frequencies (rad/s):")
for i, omega in enumerate(omega_n, 1):
    print(f"Mode {i}: {omega:.4f} rad/s")

# Time periods for each mode
T_n = 2 * np.pi / omega_n
print("\nTime periods (sec):")
for i, T in enumerate(T_n, 1):
    print(f"Mode {i}: {T:.4f} sec")


# Mode shapes (phi) before normalization
phi = eigenvectors

# Normalize each mode shape (each column) by its maximum absolute value
for n in range(num_floors):
    phi[:, n] /= np.max(np.abs(phi[:, n]))

print("\nMode shapes (columns of phi):")
for i in range(num_floors):
    print(f"Mode {i+1}: {phi[:, i]}")

# Full phi matrix view
print("\nFull Mode Shape Matrix (phi):")
print(np.round(phi, 3))  # Round for cleaner output

# Floor levels (storey numbers from 1 to num_floors)
storeys = np.arange(1, num_floors + 1)
'''
# Plot each mode shape in a separate figure with correct orientation
for i in range(num_floors):
    plt.figure(figsize=(6, 5))
    plt.plot(phi[:, i], storeys, marker='o', linestyle='-', color='b')
    plt.title(f'Mode Shape {i + 1}')
    plt.xlabel('Normalized Displacement')
    plt.ylabel('Storey Number')
    plt.grid(True)
    plt.ylim(1, num_floors)  # Bottom = Storey 1, Top = Storey 6
    plt.tight_layout()
    plt.show()
'''
# Vector of ones (for base excitation)
ones_vector = np.ones((num_floors, 1))

# Initialize containers
tau_n = []  # 1xn vector
S = np.zeros((num_floors, num_floors))  # nxn matrix (each column is s_n)

for n in range(num_floors):
    phi_n = phi[:, n].reshape(-1, 1)  # Column vector

    # Compute Gamma_n (tau_n) #Modal Participation Factor
    numerator = phi_n.T  @M@ ones_vector
    denominator = phi_n.T @ M @ phi_n
    gamma_n = (numerator / denominator).item()
    tau_n.append(gamma_n)
    EMMR = ((numerator**2)/(denominator))/((ones_vector).T@M@ones_vector)
    print("Effective Modal MAss ratio")
    print(EMMR)
    # Compute s_n and store as column in matrix S  Modal Inertial Force
    s_n = gamma_n * (M @ phi_n)
    S[:, n] = s_n.flatten()

# Display  Modal mass & tau_n vector(Moadl Participation factor)
print("Modal mass:")
print(np.round(denominator,3))
print("tau_n (modal participation factors):")
print(np.round(tau_n, 7))  # Rounded for readability

# Flip S matrix vertically so that the view is upright
S_flipped = np.flipud(S)

# Print the flipped S matrix
print("S matrix (flipped so top floor is at the top):")
print(np.round(S_flipped, 3))


# From Response Spectrum
#S_a = np.array([[1.084005], [2.166048], [2.166048], [2.166048],  [2.166048] ,[1.41264]]) # SL Value from Response spectrum
S_a = np.array([[3.294198], [6.021378], [6.021378], [6.021378],  [4.347792] ,[3.826881]]) # DL Value from Response spectrum

# Ensure S_a is a flat array for easier indexing
S_a = S_a.flatten()


#===============================================================================
# Modal displacements matrix (each column is the displacement vector u_n for mode n)
#===============================================================================
U = np.zeros((num_floors, num_floors))  # num floors X num floors

for n in range(num_floors):
    phi_n = phi[:, n].reshape(-1, 1)  # num of fllors x 1 column vector
    tau = tau_n[n]                   # scalar
    sa = S_a[n]                      # scalar
    omega_sq = eigenvalues[n]        # ω_n^2 already from eigh()

    # u_n = (phi_n * tau * sa) / ω_n²
    u_n = (phi_n * tau * sa) / omega_sq
    U[:, n] = 1000 *u_n.flatten()          # Store in U matrix

# Print all mode displacements
print("\nDisplacement vectors for each mode (columns of U):")
print(np.round(U, 16))


#===============================================================================
# DISPLACEMENT DIAGRAM VISUALIZATION
#===============================================================================
fig, axs = plt.subplots(num_floors, 1, figsize=(6, 12), sharex=True)

for n in range(num_floors):
    # Add 0 mm displacement at base
    displacements = np.insert(U[:, n], 0, 0)  # Insert 0 at base
    floors = np.arange(0, num_floors + 1)

    ax = axs[n]
    ax.plot(displacements, floors, marker='o', linewidth=2, color='b')

    ax.set_yticks(floors)
    ax.set_yticklabels(['Base'] + [f'F{i}' for i in range(1, num_floors + 1)])
    ax.set_ylim(-0.5, num_floors + 0.5)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.axvline(0, color='gray', linestyle='--', linewidth=1)

    ax.set_title(f'Mode {n+1} (T = {T_n[n]:.3f} s)', fontsize=10, loc='left')

# Global labels
fig.suptitle('Mode Shape Displacements', fontsize=14)
fig.text(0.5, 0.04, 'Displacement (mm)', ha='center', fontsize=12)
fig.text(0.04, 0.5, 'Floor Level', va='center', rotation='vertical', fontsize=12)

plt.tight_layout(rect=[0.05, 0.05, 1, 0.97])
plt.show()


#===============================================================================
#Inter storey drift
#===============================================================================

#storey drift matrix to include base (zero displacement)
storey_drift_with_base = np.zeros((num_floors, num_floors))  # 6x6

for n in range(num_floors):  # for each mode
    storey_drift_with_base[0, n] = U[0, n]  # floor 1 drift from base (assumed 0)
    for i in range(1, num_floors):
        storey_drift_with_base[i, n] = U[i, n] - U[i - 1, n]

# Print the updated storey drift matrix
print("\nStorey drift including base (rows = num of floors, columns = modes):")
print(np.round(storey_drift_with_base, 6))

# ==============================================================================
# INTER-STOREY DRIFT DIAGRAM VISUALIZATION (Floor-level values)
# ==============================================================================

fig, axs = plt.subplots(num_floors, 1, figsize=(6, 12), sharex=True)

for n in range(num_floors):
    x_values = storey_drift_with_base[:, n]              # Drift per floor
    y_values = np.arange(1, num_floors + 1)       

    ax = axs[n]
    ax.plot(x_values, y_values, marker='o', linewidth=2, color='darkorange')
    
    ax.set_yticks(y_values)
    ax.set_yticklabels([f'F{i}' for i in y_values])
    ax.set_ylim(0.5, num_floors + 0.5)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.axvline(0, color='gray', linestyle='--', linewidth=1)

    ax.set_title(f'Mode {n+1} Drift (T = {T_n[n]:.3f} s)', fontsize=10, loc='left')

# Global axis labels
fig.suptitle('Inter-storey Drift per Mode', fontsize=14)
fig.text(0.5, 0.04, 'Inter-storey Drift (mm)', ha='center', fontsize=12)
fig.text(0.04, 0.5, 'Floor Level', va='center', rotation='vertical', fontsize=12)

plt.tight_layout(rect=[0.05, 0.05, 1, 0.97])
plt.show()

#===============================================================================
# Equivalent static force matrix (correct orientation: top floor = index 0)
#===============================================================================
F_static = np.zeros((num_floors, num_floors))
for n in range(num_floors):
    F_static[:, n] = S_flipped[:, n] * S_a[n]  # already top-down from S_flipped

# Convert N to kN
F_static_kN = F_static / 1000  
# Storey shear: cumulative sum from top to bottom 
storey_shear = np.cumsum(F_static_kN, axis=0)
print("F_static in KN ,[each column is different modes from left to right]")
print(F_static_kN)  
'''
# Plot each mode's equivalent static force distribution
for n in range(num_floors):
    plt.figure(figsize=(6, 5))
    plt.barh(range(num_floors), F_static_kN[:, n], color='skyblue')
    plt.gca().invert_yaxis()  # Top floor at top
    plt.yticks(range(num_floors), [f"Floor {num_floors - i}" for i in range(num_floors)])
    plt.xlabel('Equivalent Static Force (kN)')
    plt.title(f'Mode {n + 1} Equivalent Static Force Distribution')
    plt.grid(True, axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
'''                                            
print("storey_shear in KN,[each column is different modes from left to right]")
print(storey_shear)

storey_shear_kN=storey_shear
'''
# Plotting
floor_levels = np.arange(0, num_floors + 1)
floor_labels = [f'Floor {i}' for i in range(num_floors, 0, -1)] + ['Base']

for n in range(num_floors):
    x_vals = []
    y_vals = []

    for i in range(num_floors):
        val = storey_shear_kN[i, n]
        x_vals.extend([val, val])
        y_vals.extend([i, i + 1])

    plt.figure(figsize=(8, 5))
    plt.plot(x_vals, y_vals, drawstyle='steps-pre', linewidth=2, color='tab:blue')
    plt.gca().invert_yaxis()
    plt.xlabel("Storey Shear (kN)", fontsize=12)
    plt.ylabel("Floor Level", fontsize=12)
    plt.title(f"Storey Shear Diagram - Mode {n+1} (T={T_n[n]:.3f} s)", fontsize=14)
    plt.yticks(np.arange(0, num_floors + 1), floor_labels)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()
'''
print("Base Shear (kN) for Each Mode:")
for n in range(num_floors):
    base_shear = storey_shear_kN[-1, n]  # Last row corresponds to base
    print(f"Mode {n + 1}: {base_shear:.2f} kN")
#===============================================================================
# Calculate static force per frame (divide first row by 15, others by 25)
#===============================================================================

# Create a divisor matrix
divisor_matrix = np.ones_like(storey_shear_kN) * 25  # Default divisor is 25
divisor_matrix[0, :] = 15  # First row (top floor) uses divisor 15

# Calculate force per frame
force_per_frame = storey_shear_kN=storey_shear / divisor_matrix

print("\nForce per frame (kN) [each column is different modes from left to right]:")
print(np.round(force_per_frame, 3))

#===============================================================================
# Calculate bending moments
#===============================================================================

# Half of each floor height (assuming all floors have same height)
half_h = h / 2  # 3.2m / 2 = 1.6m

# Initialize moment matrix
moment_matrix = np.zeros_like(force_per_frame)

# Calculate moments for each mode and floor
for mode in range(num_floors):
    for floor in range(num_floors):
        # Moment = Force * lever arm (half floor height)
        # For top floor (floor=0), lever arm is from top to mid-height of top floor
        # For other floors, it's the same
        moment_matrix[floor, mode] = force_per_frame[floor, mode] * half_h

print("\nBending moment (kNm) [each column is different modes from left to right]:")
print(np.round(moment_matrix, 3))

#===============================================================================
# Plotting force per frame
#===============================================================================
# Manually defined labels from Floor 5 to Base
floor_labels = ['Floor 5', 'Floor 4', 'Floor 3', 'Floor 2', 'Floor 1', 'Base']

for n in range(num_floors):
    plt.figure(figsize=(6, 5))
    plt.barh(range(num_floors), force_per_frame[:, n], color='lightgreen')
    plt.gca().invert_yaxis()  # Top floor at top
    
    # Apply labels explicitly
    plt.yticks(range(num_floors), floor_labels)
    
    plt.xlabel('Force per Frame (kN)')
    plt.title(f'Mode {n + 1} Force per Frame Distribution')
    plt.grid(True, axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()


#===============================================================================
# Model Combination by CQC and SRSS Method
#===============================================================================
# Damping ratio (needed for CQC)
zeta = 0.05  # 5% damping is typical for concrete structures

# Function to compute CQC correlation coefficient
def cqc_coeff(zeta, beta_ij):
    """Compute CQC correlation coefficient between two modes."""
    r = (8 * zeta**2 * (1 + beta_ij) * beta_ij**1.5 / 
         ((1 - beta_ij**2)**2 + 4 * zeta**2 * beta_ij * (1 + beta_ij)**2))
    return r

# Compute frequency ratios for all mode pairs
beta = np.zeros((num_floors, num_floors))
for i in range(num_floors):
    for j in range(num_floors):
        beta[i,j] = omega_n[i] / omega_n[j]  # Frequency ratio

# Initialize CQC correlation matrix
rho = np.zeros((num_floors, num_floors))
for i in range(num_floors):
    for j in range(num_floors):
        rho[i,j] = cqc_coeff(zeta, beta[i,j])

# SRSS and CQC combination of responses
def combine_modes(response_matrix, method='SRSS'):
    """Combine modal responses using SRSS or CQC."""
    num_responses = response_matrix.shape[0]
    combined = np.zeros(num_responses)
    
    if method == 'SRSS':
        for k in range(num_responses):
            combined[k] = np.sqrt(np.sum(response_matrix[k,:]**2))
    
    elif method == 'CQC':
        for k in range(num_responses):
            # Double summation for CQC
            sum_total = 0
            for i in range(num_floors):
                for j in range(num_floors):
                    sum_total += rho[i,j] * response_matrix[k,i] * response_matrix[k,j]
            combined[k] = np.sqrt(sum_total)
    
    return combined

# Apply to calculated response quantities
displacement_srss = combine_modes(U, 'SRSS')
displacement_cqc = combine_modes(U, 'CQC')

drift_srss = combine_modes(storey_drift_with_base, 'SRSS')
drift_cqc = combine_modes(storey_drift_with_base, 'CQC')

shear_srss = combine_modes(storey_shear, 'SRSS')
shear_cqc = combine_modes(storey_shear, 'CQC')

# Print results
print("\nSRSS Combined Displacements (mm):")
print(displacement_srss)
print("\nCQC Combined Displacements (mm):")
print(displacement_cqc)

print("\nSRSS Combined Drifts (mm):")
print(drift_srss)
print("\nCQC Combined Drifts (mm):")
print(drift_cqc)

print("\nSRSS Combined Shears (KN):")
print(shear_srss)
print("\nCQC Combined Shears (KN):")
print(shear_cqc)



# Create a figure with 3 subplots
plt.figure(figsize=(18, 12))

# ==============================================
# 1. Displacement Comparison
# ==============================================
plt.subplot(1, 3, 1)

# Create floor levels (including base at 0)
floor_levels = np.arange(num_floors + 1)
floor_labels = ['Base'] + [f'Floor {i}' for i in range(1, num_floors + 1)]

# Plot both methods
plt.plot(np.insert(displacement_srss, 0, 0), floor_levels, 
         'b-o', linewidth=2, markersize=8, label='SRSS')
plt.plot(np.insert(displacement_cqc, 0, 0), floor_levels, 
         'r--s', linewidth=2, markersize=8, label='CQC')


plt.xlabel('Displacement (mm)', fontsize=12)
plt.ylabel('Floor Level', fontsize=12)
plt.title('Combined Displacement', fontsize=14)
plt.legend(loc='upper right')
plt.yticks(floor_levels, floor_labels)
plt.grid(True, linestyle='--', alpha=0.7)

# ==============================================
# 2. Inter-storey Drift Comparison
# ==============================================
plt.subplot(1, 3, 2)

# Create floor levels for drift (only between floors)
drift_floor_levels = np.arange(1, num_floors + 1)
drift_floor_labels = [f'Floor {i-1}-{i}' for i in range(num_floors, 0, -1)]

# Plot both methods
plt.plot(np.flip(drift_srss), drift_floor_levels, 'b-o', linewidth=2, markersize=8, label='SRSS')
plt.plot(np.flip(drift_cqc), drift_floor_levels, 'r--s', linewidth=2, markersize=8, label='CQC')

plt.gca().invert_yaxis()
plt.xlabel('Inter-storey Drift (mm)', fontsize=12)
plt.ylabel('Floor Level', fontsize=12)
plt.title('Combined Inter-storey Drift', fontsize=14)
plt.legend(loc='upper right')
plt.yticks(drift_floor_levels, drift_floor_labels)
plt.grid(True, linestyle='--', alpha=0.7)

# ==============================================
# 3. Storey Shear Comparison (Corrected)
# ==============================================
plt.subplot(1, 3, 3)

# Function to create step plot data for shear
def create_shear_step_data(shear_values):
    x = []
    y = []
    for i in range(num_floors):
        x.extend([shear_values[i], shear_values[i]])
        y.extend([i, i + 1])
    return x, y

# Reverse shear values for bottom-up plotting
shear_srss_rev = shear_srss[::-1]
shear_cqc_rev = shear_cqc[::-1]

# Create step plot data
x_srss, y_srss = create_shear_step_data(shear_srss_rev)
x_cqc, y_cqc = create_shear_step_data(shear_cqc_rev)

# Plot SRSS and CQC lines
plt.plot(x_srss, y_srss, 'b-', linewidth=2, label='SRSS', drawstyle='steps-pre')
plt.plot(x_cqc, y_cqc, 'r--', linewidth=2, label='CQC', drawstyle='steps-pre')

# Plot markers at correct positions (top of each storey)
floor_marker_positions = np.arange(1, num_floors + 1)
plt.plot(shear_srss_rev, floor_marker_positions, 'bo', markersize=8)
plt.plot(shear_cqc_rev, floor_marker_positions, 'rs', markersize=8)

# Formatting
plt.xlabel('Storey Shear (kN)', fontsize=12)
plt.ylabel('Floor Level', fontsize=12)
plt.title('Combined Storey Shear', fontsize=14)
plt.legend(loc='upper right')

# Floor labels from base to top
floor_labels = ['Base'] + [f'Floor {i}' for i in range(1, num_floors + 1)]
plt.yticks(np.arange(0, num_floors + 1), floor_labels)

plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()


#===============================================================================
# Combined results using SRSS and CQC
#===============================================================================

# Combine force per frame using SRSS and CQC
force_per_frame_srss = combine_modes(force_per_frame, 'SRSS')
force_per_frame_cqc = combine_modes(force_per_frame, 'CQC')

# Combine moments using SRSS and CQC
moment_srss = combine_modes(moment_matrix, 'SRSS')
moment_cqc = combine_modes(moment_matrix, 'CQC')

print("\nSRSS Combined Force per Frame (kN):")
print(np.round(force_per_frame_srss, 3))
print("\nCQC Combined Force per Frame (kN):")
print(np.round(force_per_frame_cqc, 3))

print("\nSRSS Combined Bending Moment (kNm):")
print(np.round(moment_srss, 3))
print("\nCQC Combined Bending Moment (kNm):")
print(np.round(moment_cqc, 3))

#===============================================================================
# Create a comparison plot for combined results Bending moment
#===============================================================================
floor_levels = np.arange(num_floors)

# Manually defined floor labels from top to base
floor_labels = ['Floor 5', 'Floor 4', 'Floor 3', 'Floor 2', 'Floor 1', 'Base']

plt.figure(figsize=(6, 5))
plt.plot(force_per_frame_srss, floor_levels, 'b-o', linewidth=2, markersize=8, label='SRSS')
plt.plot(force_per_frame_cqc, floor_levels, 'r--s', linewidth=2, markersize=8, label='CQC')
plt.gca().invert_yaxis()
plt.xlabel('Force per Frame (kN)', fontsize=12)
plt.ylabel('Floor Level', fontsize=12)
plt.title('Combined Force per Frame', fontsize=14)
plt.legend(loc='upper right')

# Apply the custom floor labels
plt.yticks(floor_levels, floor_labels)

plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
