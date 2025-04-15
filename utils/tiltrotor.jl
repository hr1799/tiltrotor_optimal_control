T = Diagonal([1; -ones(3)])
H = [zeros(1,3); I]

function dcm_from_mrp(p)
    p1,p2,p3 = p
    den = (p1^2 + p2^2 + p3^2 + 1)^2
    a = (4*p1^2 + 4*p2^2 + 4*p3^2 - 4)
    return [
    (-((8*p2^2+8*p3^2)/den-1)*den)   (8*p1*p2 + p3*a)     (8*p1*p3 - p2*a);
    (8*p1*p2 - p3*a) (-((8*p1^2 + 8*p3^2)/den - 1)*den)   (8*p2*p3 + p1*a);
    (8*p1*p3 + p2*a)  (8*p2*p3 - p1*a)  (-((8*p1^2 + 8*p2^2)/den - 1)*den)
    ]/den
end

function skew(ω::Vector)
    return [0    -ω[3]  ω[2];
            ω[3]  0    -ω[1];
           -ω[2]  ω[1]  0]
end

# quaternion utilities
function L_op(q)
    s = q[1]
    v = q[2:4]
    L = [s    -v';
         v  s*I+skew(v)]
    return L
end

function qtoQ(q)
    return H'*T*L_op(q)*T*L_op(q)*H
end

function G(q)
    return G = L_op(q)*H
end

function rptoq(ϕ)
    return (1/sqrt(1+ϕ'*ϕ))*[1; ϕ]
end

function qtorp(q)
    return q[2:4]/q[1]
end

# Coefficient of lift and drag for the wing (DAE 51 airfoil)
function C_l(α)
    return -7e-6*α^3 - 0.0013*α^2 + 0.0563*α + 0.4683
end

function C_d(α)
    return -4e-6*α^3 - 0.0005*α^2 + 0.0008*α + 0.0269
end

function tiltrotor_dynamics_mrp(model::NamedTuple,x,u)
    # tiltrotor dynamics with an MRP for attitude
    # and velocity in the world frame (not body frame)
    
    r = x[1:3]     # position in world frame 
    v = x[4:6]     # velocity in body frame 
    p = x[7:9]     # n_p_b (MRP) attitude 
    ω = x[10:12]   # angular velocity 
    δ_l = x[13]    # left rotor angle
    δ_r = x[14]    # right rotor angle

    Q = dcm_from_mrp(p) # Rotation matrix from body to world frame

    mass = model.mass
    J = model.J
    gravity = model.gravity
    L = model.L
    kf = model.kf
    km = model.km
    kδ = model.kδ

    w1 = u[1]
    w2 = u[2]

    F_l = max(0,kf*w1)
    F_r = max(0,kf*w2)
    M_l = km*w1
    M_r = -km*w2

    δ_l_dot = kδ*u[3]
    δ_r_dot = kδ*u[4]
    
    # Force due to propellers
    F_thrusters = [F_r*sin(δ_r) + F_l*sin(δ_l); 0; F_r*cos(δ_r) + F_l*cos(δ_l)]

    # Force due to wings
    α = -atan(v[3], v[1]) #Angle of attack
    q = 0.5*model.ρ*(v[1]^2 + v[3]^2)
    L_wing = [0; 0; q*model.S_wing*C_l(α)] #Lift force
    D_wing = [-q*model.S_wing*C_d(α); 0; 0] #Drag force

    F = F_thrusters + L_wing + D_wing

    f = Q'*mass*gravity + F # forces in world frame

    # Moments due to propellers
    M_thrusters = [M_r*sin(δ_r) + M_l*sin(δ_l); 0; M_r*cos(δ_r) + M_l*cos(δ_l)]
    
    τ = [L*(-F_r*cos(δ_r) + F_l*cos(δ_l)), 0, (F_r*sin(δ_r) - F_l*sin(δ_l))] + M_thrusters #total rotor torque in body frame

    # this is xdot 
    return [
        Q*v
        f/mass - skew(ω)*v
        ((1+norm(p)^2)/4)*(I + 2*(skew(p)^2 + skew(p))/(1+norm(p)^2))*ω
        J\(τ - cross(ω,J*ω))
        δ_l_dot
        δ_r_dot
    ]
end

function tiltrotor_dynamics_quat(model::NamedTuple,x,u)
    # tiltrotor dynamics with a quaternion for attitude
    
    r = x[1:3]                  # position in world frame 
    v = x[4:6]                  # velocity in body frame 
    quat = x[7:10]/norm(x[7:10])   # normalize q just to be careful
    ω = x[11:13]                # angular velocity 
    δ_l = x[14]                 # left rotor angle
    δ_r = x[15]                 # right rotor angle

    Q = qtoQ(quat) # Rotation matrix from body to world frame

    mass = model.mass
    J = model.J
    gravity = model.gravity
    L = model.L
    kf = model.kf
    km = model.km
    kδ = model.kδ

    w1 = u[1]
    w2 = u[2]

    F_l = max(0,kf*w1)
    F_r = max(0,kf*w2)
    M_l = km*w1
    M_r = -km*w2

    δ_l_dot = kδ*u[3]
    δ_r_dot = kδ*u[4]

    # Force due to propellers
    F_thrusters = [F_r*sin(δ_r) + F_l*sin(δ_l); 0; F_r*cos(δ_r) + F_l*cos(δ_l)]

    # Force due to wings
    α = -atan(v[3], v[1]) #Angle of attack
    q = 0.5*model.ρ*(v[1]^2 + v[3]^2)
    L_wing = [0; 0; q*model.S_wing*C_l(α)] #Lift force
    D_wing = [-q*model.S_wing*C_d(α); 0; 0] #Drag force

    F = F_thrusters + L_wing + D_wing

    f = Q'*mass*gravity + F # forces in world frame

    # Moments due to propellers
    M_thrusters = [M_r*sin(δ_r) + M_l*sin(δ_l); 0; M_r*cos(δ_r) + M_l*cos(δ_l)]
    
    τ = [L*(-F_r*cos(δ_r) + F_l*cos(δ_l)), 0, (F_r*sin(δ_r) - F_l*sin(δ_l))] + M_thrusters #total rotor torque in body frame

    # this is xdot 
    return [
        Q*v
        f/mass - skew(ω)*v
        0.5*L_op(quat)*H*ω
        J\(τ - cross(ω,J*ω)) # J\(τ - skew(ω)*J*ω) 
        δ_l_dot
        δ_r_dot
    ]
end

function rk4(model,ode,x,u,dt)
    # rk4 
    k1 = dt*ode(model,x, u)
    k2 = dt*ode(model,x + k1/2, u)
    k3 = dt*ode(model,x + k2/2, u)
    k4 = dt*ode(model,x + k3, u)
    return x + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
end

function hermite_simpson(model::NamedTuple, ode, x1::Vector, x2::Vector, u, dt::Real)::Vector
    f1 = ode(model, x1, u)
    f2 = ode(model, x2, u)
    xm = 0.5*(x1 + x2) + (dt/8.0)*(f1 - f2)
    ẋm = (-3/(2.0*dt))*(x1 - x2) - 0.25*(f1 + f2)
    fm = ode(model, xm, u)
    return fm - ẋm
end

function animate_tiltrotor_mrp(X, dt)
    vis = mc.Visualizer()
    urdf = joinpath(@__DIR__,"../urdf/tiltrotor.urdf")
    robot = parse_urdf(urdf)
    remove_fixed_tree_joints!(robot)

    mvis = MechanismVisualizer(robot, URDFVisuals(urdf), vis)
    # settransform!(vis, Translation(0., 0, 3))
    # set_configuration!(mvis, [1.5, -1.5])
    anim = mc.Animation(floor(Int,1/dt))
    for k = 1:length(X)
        mc.atframe(anim, k) do
            p = X[k][7:9]
            mc.settransform!(vis, Translation(X[k][1], X[k][2], X[k][3]) ∘ LinearMap((dcm_from_mrp(p))))
            set_configuration!(mvis, [X[k][13], X[k][14]])
        end
    end
    mc.setanimation!(vis, anim)
    return mc.render(vis)
end

function animate_tiltrotor_quat(X, dt)
    vis = mc.Visualizer()
    urdf = joinpath(@__DIR__,"../urdf/tiltrotor.urdf")
    robot = parse_urdf(urdf)
    remove_fixed_tree_joints!(robot)

    mvis = MechanismVisualizer(robot, URDFVisuals(urdf), vis)
    anim = mc.Animation(floor(Int,1/dt))
    for k = 1:length(X)
        mc.atframe(anim, k) do
            q = QuatRotation(X[k][7], X[k][8], X[k][9], X[k][10])
            mc.settransform!(vis, Translation(X[k][1], X[k][2], X[k][3]) ∘ LinearMap(qtoQ(q)))
            set_configuration!(mvis, [X[k][14], X[k][15]])
        end
    end
    mc.setanimation!(vis, anim)
    return mc.render(vis)
end
