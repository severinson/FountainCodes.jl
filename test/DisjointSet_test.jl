using FountainCodes, Test

function test_1()
    st = IntDisjointSetsTracked(10)
    union!(st, 1, 2)
    if cvertices(st) != 2
        error("cvertices = $(cvertices(st)) != 2")
    end
    if cedges(st) != 1
        error("cedges = $(cedges(st)) != 1")
    end
    if !(croot(st) in [1, 2])
        error("croot = $(croot(st)) not in [1, 2]")
    end
    return true
end
@test test_1()

function test_2()
    st = IntDisjointSetsTracked(10)
    union!(st, 1, 2)
    union!(st, 2, 3)
    union!(st, 3, 4)
    union!(st, 4, 5)
    union!(st, 5, 6)
    union!(st, 6, 7)
    cv = 7
    if cvertices(st) != cv
        error("cvertices = $(cvertices(st)) != $cv")
    end
    ce = 6
    if cedges(st) != ce
        error("cedges = $(cedges(st)) != $ce")
    end
    union!(st, 1, 7)
    cv = 7
    if cvertices(st) != cv
        error("cvertices = $(cvertices(st)) != $cv")
    end
    ce = 7
    if cedges(st) != ce
        error("cedges = $(cedges(st)) != $ce")
    end
    if !(croot(st) in [1, 2, 3, 4, 5, 6, 7])
        error("croot = $(croot(st)) not in [1, 2]")
    end
    return true
end
@test test_2()

function test_3()
    st = IntDisjointSetsTracked(20)
    union!(st, 1, 2)
    union!(st, 2, 3)
    union!(st, 3, 4)
    union!(st, 4, 5)
    union!(st, 5, 1)
    cv = 5
    if cvertices(st) != cv
        error("cvertices = $(cvertices(st)) != $cv")
    end
    ce = 5
    if cedges(st) != ce
        error("cedges = $(cedges(st)) != $ce")
    end

    union!(st, 6, 7)
    union!(st, 7, 8)
    union!(st, 8, 9)
    union!(st, 10, 11)
    union!(st, 11, 12)
    union!(st, 12, 6)
    cv = 7
    if cvertices(st) != cv
        error("cvertices = $(cvertices(st)) != $cv")
    end
    ce = 6
    if cedges(st) != ce
        error("cedges = $(cedges(st)) != $ce")
    end

    union!(st, 4, 9)
    cv = 12
    if cvertices(st) != cv
        error("cvertices = $(cvertices(st)) != $cv")
    end
    ce = 12
    if cedges(st) != ce
        error("cedges = $(cedges(st)) != $ce")
    end
    if !(croot(st) in collect(1:12))
        error("croot = $(croot(st)) not in $(collect(1:12))")
    end
    return true
end
@test test_3()

function test_4()
    st = IntDisjointSetsTracked(10)
    union!(st, 1, 2)
    union!(st, 1, 3)
    cv = 3
    if cvertices(st) != cv
        error("cvertices = $(cvertices(st)) != $cv")
    end
    ce = 2
    if cedges(st) != ce
        error("cedges = $(cedges(st)) != $ce")
    end
    if !(croot(st) in [1, 2, 3])
        error("croot = $(croot(st)) not in [1, 2, 3]")
    end

    reset!(st)
    union!(st, 1, 2)
    union!(st, 1, 3)
    cv = 3
    if cvertices(st) != cv
        error("cvertices = $(cvertices(st)) != $cv")
    end
    ce = 2
    if cedges(st) != ce
        error("cedges = $(cedges(st)) != $ce")
    end
    if !(croot(st) in [1, 2, 3])
        error("croot = $(croot(st)) not in [1, 2, 3]")
    end
    return true
end
@test test_4()

function test_5()
    st = IntDisjointSetsTracked(10)
    v = 1
    if vertices(st, 1) != v
        error("vertices = $(vertices(st, 1)) != $v")
    end
    e = 0
    if edges(st, 1) != e
        error("edges = $(edges(st, 1)) != $e")
    end

    union!(st, 1, 2)
    union!(st, 1, 3)
    v = 3
    if vertices(st, 3) != v
        error("vertices = $(vertices(st, 1)) != $v")
    end
    e = 2
    if edges(st, 3) != e
        error("edges = $(edges(st, 1)) != $e")
    end

    # cv = 3
    # if cvertices(st) != cv
    #     error("cvertices = $(cvertices(st)) != $cv")
    # end
    # ce = 2
    # if cedges(st) != ce
    #     error("cedges = $(cedges(st)) != $ce")
    # end
    # if !(croot(st) in [1, 2, 3])
    #     error("croot = $(croot(st)) not in [1, 2, 3]")
    # end

    # reset!(st)
    # union!(st, 1, 2)
    # union!(st, 1, 3)
    # cv = 3
    # if cvertices(st) != cv
    #     error("cvertices = $(cvertices(st)) != $cv")
    # end
    # ce = 2
    # if cedges(st) != ce
    #     error("cedges = $(cedges(st)) != $ce")
    # end
    # if !(croot(st) in [1, 2, 3])
    #     error("croot = $(croot(st)) not in [1, 2, 3]")
    # end
    return true
end
@test test_5()
