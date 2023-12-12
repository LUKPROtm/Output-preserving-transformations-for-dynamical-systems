#### GENERIC FUNCTION TO FIND AUTOMATICALLY SYMBOLIC SOLUTIONS (FOR SIMPLE MODELS)#######################################################

using SymPy
using AbstractAlgebra
using HomotopyContinuation
using StructuralIdentifiability



function find_transformations(model, names_map, with_states = false)
    id_funct = find_identifiable_functions(model , with_states = with_states)
    for var in names_map
        @vars names_map[1], names_map[2]
    end
    @vars u
    system = create_system_of_equations(id_funct, names_map) 
    unknown = [var[2] for var in names_map]
    SymPy.solve(system, unknown)
end


#### CREATE SYMBOLIC SYSTEM OF EQUATIONS FROM IDENTIFIABLE FUNCTIONS#####################################################################

function create_system_of_equations(id_funct, names_map)
    id_funct_num_sympy = []
    id_funct_den_sympy = []
    for el in id_funct
        num = numerator(el)
        den = denominator(el)
        push!(id_funct_num_sympy, sympify(string(num)))
        push!(id_funct_den_sympy, sympify(string(den)))
    end
    
    id_funct_num_sympy_transf = copy(id_funct_num_sympy)
    id_funct_den_sympy_transf = copy(id_funct_den_sympy)
    
    for j in 1:length(id_funct_num_sympy_transf)
        id_funct_num_sympy_transf[j] = id_funct_num_sympy_transf[j].subs(names_map)
    end
    for j in 1:length(id_funct_den_sympy_transf)
            id_funct_den_sympy_transf[j] = id_funct_den_sympy_transf[j].subs(names_map)
    end
    
    lhs = id_funct_num_sympy .* id_funct_den_sympy_transf
    rhs = id_funct_num_sympy_transf .* id_funct_den_sympy
    system = lhs - rhs
    return system
    
end


#### CREATE SYMBOLIC SYSTEM OF EQUATIONS FROM IDENTIFIABLE FUNCTIONS ALGEBRA #####################################################################

#to be done

function create_system_of_equations_nemo(id_funct, names_map)
    id_funct_num_sympy = []
    id_funct_den_sympy = []
    for el in id_funct
        num = numerator(el)
        den = denominator(el)
        push!(id_funct_num_sympy, sympify(string(num)))
        push!(id_funct_den_sympy, sympify(string(den)))
    end
    
    id_funct_num_sympy_transf = copy(id_funct_num_sympy)
    id_funct_den_sympy_transf = copy(id_funct_den_sympy)
    
    for j in 1:length(id_funct_num_sympy_transf)
        id_funct_num_sympy_transf[j] = id_funct_num_sympy_transf[j].subs(names_map)
    end
    for j in 1:length(id_funct_den_sympy_transf)
            id_funct_den_sympy_transf[j] = id_funct_den_sympy_transf[j].subs(names_map)
    end
    
    lhs = id_funct_num_sympy .* id_funct_den_sympy_transf
    rhs = id_funct_num_sympy_transf .* id_funct_den_sympy
    system = lhs - rhs
    return system
    
end
############################### UNUSED #################################
function create_system_of_equations_num_sympy(id_funct, names, values)
    id_funct_num_sympy = []
    id_funct_den_sympy = []
    for el in id_funct
        num = numerator(el)
        den = denominator(el)
        push!(id_funct_num_sympy, sympify(string(num)))
        push!(id_funct_den_sympy, sympify(string(den)))
    end
    
    id_funct_num_sympy_transf = copy(id_funct_num_sympy)
    id_funct_den_sympy_transf = copy(id_funct_den_sympy)
    
    map = [(names[i], values[i]) for i in 1:length(names)]
    
    for j in 1:length(id_funct_num_sympy_transf)
        id_funct_num_sympy_transf[j] = id_funct_num_sympy_transf[j].subs(map)
    end
    for j in 1:length(id_funct_den_sympy_transf)
            id_funct_den_sympy_transf[j] = id_funct_den_sympy_transf[j].subs(map)
    end
    
    lhs = id_funct_num_sympy .* id_funct_den_sympy_transf
    rhs = id_funct_num_sympy_transf .* id_funct_den_sympy
    system = lhs - rhs
    return system
    
end

######################### CREATE NUMERICAL SYSTEM OF EQUATIONS TO BE USED WITH HomotopyContinuation###########################################

function create_system_of_equations_num_hc(id_funct, names, initial_conditions)
    id_funct_num_hc = Array{Expression}(undef, 0)
    id_funct_den_hc = Array{Expression}(undef, 0)
    for el in id_funct
        num = numerator(el)
        den = denominator(el)
        push!(id_funct_num_hc, nemo2hc(Meta.parse(string(num))))
        push!(id_funct_den_hc, nemo2hc(Meta.parse(string(den))))
    end
	 
    # Substitute variable values into the polynomials
    id_funct_num_hc_transf = [HomotopyContinuation.subs(p, names => initial_conditions) for p in id_funct_num_hc]
    id_funct_den_hc_transf = [HomotopyContinuation.subs(p, names => initial_conditions) for p in id_funct_den_hc]
	
    # Create systems of equations
    lhs = id_funct_num_hc .* id_funct_den_hc_transf
    rhs = id_funct_num_hc_transf .* id_funct_den_hc
    system = lhs .- rhs

    return system
    
end


function find_numerical_transformations_sympy(model, names, with_states = false)
    id_funct = find_identifiable_functions(model , with_states = with_states)
    for var in names
        @vars eval(names_map[1])
    end
    @vars u
    values = rand(Float64, length(names))
    system = create_system_of_equations_num_sympy(id_funct,names, values) 
    SymPy.solve(system)
end


function find_numerical_transformations_hc(model, names, initial_conditions, with_states = false)
    id_funct = find_identifiable_functions(model , with_states = with_states)
    for var in names
        @var names_map[1]
    end
    @var u
	varnames= Array{Variable}(undef, 0)
	for name in names
		push!(varnames, Variable(Meta.parse(name)))
	end
	#push!(varnames, u)
 
    system = create_system_of_equations_num_hc(id_funct,varnames, initial_conditions) 
    return results(HomotopyContinuation.solve(system), only_real = true)
end

function nemo2hc(expr_tree::Union{Expr,Symbol, Int64})
    #traverse expr_tree
	if typeof(expr_tree) == Int64
        return Expression(expr_tree)
    end
    if typeof(expr_tree) == Symbol
        return HomotopyContinuation.variables(expr_tree)[1]
    end
    if typeof(expr_tree) == Expr
        if expr_tree.head == :call
            if expr_tree.args[1] in [:+, :-, :*, :/, :^]
                return reduce(eval(expr_tree.args[1]), map(nemo2hc, expr_tree.args[2:end]))
            end
        end
    end
end
