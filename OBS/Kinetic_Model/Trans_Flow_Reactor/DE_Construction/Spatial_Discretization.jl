function first_order_upwind(dx::Vector,upwind='L')
    # this version for variable spacing
    if upwind=='L'
        Mat=LinearAlgebra.diagm(0=>1 ./ dx, -1=> -1 ./ dx[2:end])
    else
        Mat=Mat=LinearAlgebra.diagm(0=>-1 ./ dx, 1=>1 ./ dx[1:end-1])
    end
    return Mat
end

function first_order_upwind(dx::AbstractFloat, N::Int, upwind='L')
    # this version is for uniform spacing
    return @inline first_order_upwind(dx*ones(Float64, N), upwind)
end

function central_diff(dx::AbstractFloat)
    # WARNING: I got lazy and this only works for uniform spacing
    # I try and enforce this by only accepting a float
    # switches to 1st order on boundaries
    Mat=SparseArrays.spdiagm(0=>[-1/dx;zeros(NvNc[2]-2);1/dx],
        -1=>[-ones(NvNc[2]-2) ./ (2*dx);-1/dx],
        1=>[1/dx;ones(NvNc[2]-2) ./ (2*dx)])
    return Mat
end