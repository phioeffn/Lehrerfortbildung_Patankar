using LinearAlgebra
using DifferentialEquations
using Plots



function patankar(prod_dest, tspan, u0)

    dim = length(u0)             # Dimension of the problem
    Nt = length(tspan)           # Length of time span
    U = zeros(dim, Nt)           # Solution vector
    p = zeros(dim, dim)          # Temporary production matrix
    d = zeros(dim, dim)          # Temporary destruction matrix
    U[:, 1] = u0                 # Set initial condition (Julia is 1-based indexing)

    for it in 2:Nt               # Loop over timesteps (1-based indexing in Julia)
        dt = tspan[it] - tspan[it-1]   # Time step size
        p, d = prod_dest(U[:, it-1])   # Compute production and destruction at previous timestep
        
        for i in 1:dim                 # Loop over each dimension
            lhs = 1.0                  # Initialize lhs
            rhs = U[i, it-1]           # Initialize rhs with the previous value
            for j in 1:dim             # Loop over each element in matrices p and d
                lhs += dt * d[i,j] / U[i, it-1]  # Update lhs
                rhs += dt * p[i,j]               # Update rhs
            end
           
            U[i, it] = rhs / lhs        # Solve the final system for U
         end
    end
    return tspan, U
end


function mPEuler(prod_dest, tspan, u0)
    """
    Input:
    - prod_dest: function that returns the matrices p[i,j](c) and d[i,j](c)
    - tspan: time vector
    - u0: initial condition
    """
    
    dim = length(u0)             # Dimension of the problem
    Nt = length(tspan)           # Length of time span
    U = zeros(dim, Nt)           # Solution vector
    p = zeros(dim, dim)          # Temporary production matrix
    d = zeros(dim, dim)          # Temporary destruction matrix
    U[:, 1] = u0                 # Set initial condition (Julia is 1-based indexing)

    for it in 2:Nt               # Loop over timesteps (1-based indexing in Julia)
        dt = tspan[it] - tspan[it-1]   # Time step size
        p, d = prod_dest(U[:, it-1])   # Compute production and destruction at previous timestep
        
        MM = Matrix{Float64}(I, dim, dim)  # Initialize the mass matrix (identity matrix)

        for i in 1:dim                 # Loop over each dimension
            for j in 1:dim             # Loop over each element in matrices p and d
                MM[i,j] -= dt * p[i,j] / U[j, it-1]  # Update off-diagonal entries of MM
                MM[i,i] += dt * d[i,j] / U[i, it-1]  # Update diagonal entries of MM
            end
        end

        # Solve the system: MM * U[:,it] = U[:,it-1]
        U[:, it] = MM \ U[:, it-1]      # Use Julia's backslash operator to solve linear system
    end

    return tspan, U
end