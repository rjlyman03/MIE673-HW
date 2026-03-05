import numpy as np
import matplotlib.pyplot as plt


# Gaussian Markov Chain

def generate_ar1(alpha, N):
    """
    Generates a time series of length N using the AR(1) / Gaussian Markov Chain:
        X_{i+1} = alpha * X_i + Z_i
    where Z_i ~ N(0, 1)
    and   X_1  ~ N(0, 1/(1 - alpha^2))
    """
    X = np.zeros(N)  # pre-allocate the array

    # --- X_1: draw from the stationary distribution N(0, 1/(1-alpha^2))
    X[0] = np.random.normal(loc=0, scale=np.sqrt(1 / (1 - alpha**2)))

    # --- X_2 ... X_N: apply the recurrence relation
    for i in range(1, N):
        Z = np.random.normal(0, 1)      # independent noise term
        X[i] = alpha * X[i-1] + Z      # AR(1) update rule

    return X



# alpha values and series length

alpha_values = [0.2, 0.5, 0.9, 0.99]   
N = 1000                                 

np.random.seed(42)  



# plot each time series as individual image

for alpha in alpha_values:
    X = generate_ar1(alpha, N)

    #  variance from the formula: σ² = 1 / (1 - α²)
    theoretical_var = 1 / (1 - alpha**2)
    empirical_var   = np.var(X)

    
    fig, ax = plt.subplots(figsize=(12, 4))

    ax.plot(X, linewidth=0.8, color='steelblue')
    ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.set_title(
        f"Gaussian Markov Chain  |  α = {alpha}\n" ,
        #f"Theoretical σ² = {theoretical_var:.3f}  |  Empirical σ² = {empirical_var:.3f}",
        fontsize=12, fontweight='bold'
    )
    ax.set_xlabel("Time step i", fontsize=11)
    ax.set_ylabel("X_i", fontsize=11)

    plt.tight_layout()

    # Save each plot as its own file
    filename = f"ar1_alpha_{str(alpha).replace('.', '')}.png"
    plt.savefig(filename, dpi=150)
    plt.close()  

    print(f"Saved: {filename}")


# (b) autocorrelation, exponential decay

# define autocorrelation function
def autocorrelation(X, max_lag):
    """
    Computes the autocorrelation of series X for lags 0, 1, 2, ..., max_lag.
    
    Autocorrelation at lag k = Corr(X_i, X_{i+k})
    = E[(X_i - mean)(X_{i+k} - mean)] / variance
    """
    X_centered = X - np.mean(X)        # remove the mean
    var = np.var(X)                    # variance of the series
    N = len(X)

    acf = []
    for k in range(max_lag + 1):
        # multiply each value by the value k steps ahead, then average
        cov_k = np.mean(X_centered[:N-k] * X_centered[k:])
        acf.append(cov_k / var)        # normalize → correlation (between -1 and 1)

    return np.array(acf)


# set parameters and compute autocorrelation
alpha_values = [0.2, 0.5, 0.9, 0.99]
N       = 5000          # longer series goes, smoother autocorrelation estimate
max_lag = 50            # how many lags to look at

np.random.seed(42)

# plot autocorrelation vs exponential decay for each alpha
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.flatten()   # makes it easy to loop over

fig.suptitle("Autocorrelation vs Exponential Decay  (ACF ≈ αᵏ)", fontsize=14, fontweight='bold')

lags = np.arange(0, max_lag + 1)   # [0, 1, 2, ..., max_lag]

for ax, alpha in zip(axes, alpha_values):
    # --- empirical autocorrelation from the generated data
    X   = generate_ar1(alpha, N)          # reuse from part (a)
    acf = autocorrelation(X, max_lag)

    #  exponential decay: R(k) = alpha^k
    # theoretical autocorrelation for an AR(1) process/ Gaussian Markov Chain
    exp_decay = alpha ** lags

    # plotting both
    ax.plot(lags, acf,       color='steelblue', linewidth=1.5, label='Empirical ACF')
    ax.plot(lags, exp_decay, color='tomato',    linewidth=2,
            linestyle='--',  label=f'Analytical: α^k = {alpha}^k')

    ax.set_title(f"α = {alpha}", fontsize=11)
    ax.set_xlabel("Lag  k")
    ax.set_ylabel("Autocorrelation  R(k)")
    ax.legend(fontsize=8)
    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')
    ax.set_ylim(-0.2, 1.05)

plt.tight_layout()
plt.savefig("ar1_autocorrelation.png", dpi=150)
plt.show()
