
## Linear Example 

function prod_dest(c)
    d= length(c)
    p=zeros(d,d)
    d=zeros(d,d)
    p[1,2]=param_a*c[2]
    d[2,1]=param_a*c[2]
    p[2,1]=param_b*c[1]
    d[1,2]=param_b*c[1]
    return p, d
end




# Define the SIR flux (derivatives) function
function linear_system2_flux!(du, u, p, t)
    beta,gamma = p  # Extract parameters (infection rate and recovery rate)
    du[1] = -5 * u[1] + u[2]  # dS/dt: rate of change of susceptible individuals
    du[2] =  5 * u[1] - u[2]  # dI/dt: rate of change of infected individuals
end

function linear_system2_jacobian!(J, u, p, t)
    # Fill the Jacobian matrix J with the partial derivatives
    J[1, 1] = -5   # ∂(dS/dt)/∂u[1]
    J[1, 2] =  1   # ∂(dS/dt)/∂u[2]
    J[2, 1] =  5   # ∂(dI/dt)/∂u[1]
    J[2, 2] = -1   # ∂(dI/dt)/∂u[2]
end






#### SIR-Modell 

function SIR_production_destruction(u, t=0, beta=1.0, gamma=0.1)
    dim = length(u)
    p = zeros(dim, dim)  # Production matrix
    d = zeros(dim, dim)  # Destruction matrix
    N = sum(u)           # Total population (S + I + R)

    # Production and destruction terms
    p[2, 1] = beta * u[1] * u[2] / N  # Production term: S -> I
    d[1, 2] = p[2, 1]                 # Destruction term: I -> S

    p[3, 2] = gamma * u[2]            # Production term: I -> R
    d[2, 3] = p[3, 2]                 # Destruction term: R -> I

    return p, d
end




# Define the SIR flux (derivatives) function
function SIR_flux!(du, u, p, t)
    beta,gamma = p  # Extract parameters (infection rate and recovery rate)
    N = sum(u)       # Total population (S + I + R)
    
    du[1] = -beta * u[1] * u[2] / N  # dS/dt: rate of change of susceptible individuals
    du[2] = beta * u[1] * u[2] / N - gamma * u[2]  # dI/dt: rate of change of infected individuals
    du[3] = gamma * u[2]  # dR/dt: rate of change of recovered individuals
end



#### Roberts-Test 


function rob_prod_dest(c)
    d= length(c)
    p=zeros(d,d)
    d=zeros(d,d)
    p[1,2]=10^4*c[3]^1.5 #"c[2]"
    d[2,1]=10^4*c[3]^1.5 #"c[2]"
    p[2,1]=0.04*c[1]
    d[1,2]=0.04*c[1]
    p[3,2]=3*10^7*c[2]^1.5 #"^2"
    d[2,3]=3*10^7*c[2]^1.5 #"^2"
    return p, d
 end


 function Robertson_production_destruction(u; alpha=1e4, beta=0.04, gamma=3e7)
        # u is the vector of current state variables
        # alpha, beta, and gamma are the parameters
    
        dim = length(u)
        p = zeros(dim, dim)  # Production matrix
        d = zeros(dim, dim)  # Destruction matrix
    
        # Populate the production and destruction matrices
        p[1, 2] = alpha * u[2] * u[3]
        d[2, 1] = p[1, 2]
        
        p[2, 1] = beta * u[1]
        d[1, 2] = p[2, 1]
        
        p[3, 2] = gamma * u[2]^2
        d[2, 3] = p[3, 2]
    
        return p, d
    end


   # Define the Robertson ODE system
   function Robertson!(du, u, p, t)
    y1, y2, y3 = u
    du[1] = -0.04 * y1 + 1e4 * y2 * y3
    du[2] = 0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2^2
    du[3] = 3e7 * y2^2
end



### 




