grouped_xy(x::AbstractVector, y::AbstractMatrix) = x, y
grouped_xy(y::AbstractMatrix) = 1:size(y,1), y

@recipe function f(cells::Array{CINulator.cancercellCN, 1})

    f = copynumberfrequency(cells)
    x, y = grouped_xy(1:cells[1].chromosomes.N[1], f')
    mylabels = map(x -> "CN = $x", 0:length(f[:, 1]))
    x = map(x -> "$x", x)


    nr, nc = size(y)
    isstack = pop!(plotattributes, :bar_position, :dodge) == :stack

    xnums = if eltype(x) <: Number
        bar_width --> (0.8 * mean(diff(x)))
        x
    else
        bar_width --> 0.8
        ux = unique(x)
        xnums = (1:length(ux)) .- 0.5
        xticks --> (xnums, ux)
        lab --> mylabels
        xlab --> "Chromosome"
        palette --> :inferno
        xnums
    end
    @assert length(xnums) == nr


    # compute the x centers.  for dodge, make a matrix for each column
    x = if isstack
        x
    else
        bws = plotattributes[:bar_width] / nc
        bar_width := bws
        xmat = zeros(nr,nc)
        for r=1:nr
            bw = _cycle(bws, r)
            farleft = xnums[r] - 0.5 * (bw * nc)
            for c=1:nc
                xmat[r,c] = farleft + 0.5bw + (c-1)*bw
            end
        end
        xmat
    end

    # compute fillrange
    fillrange := if isstack
        # shift y/fillrange up
        y = copy(y)
        y[.!isfinite.(y)] .= 0
        fr = zeros(nr, nc)
        for c=2:nc
            for r=1:nr
                fr[r,c] = y[r,c-1]
                y[r,c] += fr[r,c]
            end
        end
        fr
    else
        get(plotattributes, :fillrange, nothing)
    end

    seriestype := :bar
    x, y
end
