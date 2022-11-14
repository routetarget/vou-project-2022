using JuMP, HiGHS

## Data
begin
    #TODO dalsi makra na rozmery matic tady
    macro NUM_SCENARIOS()
        RETURN :(10)
    end
    macro NUM_ITEMS()
        RETURN :(15)
    end
    macro NUM_CENTERS()
        RETURN :(60)
    end
    macro NUM_COV_LEVELS()
        RETURN :(5)
    end

    #TODO generace dat v maticich, zatim hotove jen velikosti matic
    #TODO prejmenovat promenne na citelne nazvy
    P = zeros((1,@NUM_SCENARIOS)); # Matice pravdepodobnosti udalosi s
    D = zeros((@NUM_SCENARIOS, @NUM_ITEMS)); # matice vsech ocekavanych poptavek po itemu k v udalosi S
    Cap = zeros((1, @NUM_CENTERS)); # kapacita distribucnich center (volume - objem)
    B_0 = 1200000; # rozpocet pred katastrofou
    B_1 = 2*B_0; # rozpocet po katastrofe
    F = zeros((1, @NUM_CENTERS)); # Budget pro zavedeni konkretniho centra
    G = zeros((@NUM_CENTERS, @NUM_ITEMS)); # cena itamu typu k v distribucnim centru j ($/jednotku)
    C = zeros((@NUM_CENTERS, @NUM_ITEMS, @NUM_SCENARIOS)); # cena dopravy itemu k z centra j na lokaci katastrofy s ($/jednotku)
    T = zeros((@NUM_CENTERS, @NUM_ITEMS, @NUM_SCENARIOS)); # Cas naplneni dodavky itamu k z centra j na katastrofu j (hodiny)
    W = zeros((1, @NUM_ITEMS)); # Zavasnost pro item k, sum_k(W_k) = 1, W_k >= 0
    L = zeros((1, @NUM_ITEMS)); #level pokryti pro item k
    A = zeros((@NUM_ITEMS)); # Hodnota levelu pokryti, mezi 1 a 0, serazeno, kazda odpovida jednomu levelu pokryti
    N_s = # Mnozina dist. center ktere jsou schopny splnit level pokryti lk pro item typu k v udalosi s
    LR = # spodni limit jednoho levelu pokryti lk
    UR = # horni limit jednoho levelu pokryti lk



end

## model
begin
    model = Model(HiGHS.Optimizer)
    # Decicision variables
    @variable(model, f) # Podil itemu typu k, naplnen distribucnim centrem k, v udalosti s
    @variable(model, Q) # Kusu itemu typu k uskladneneho v distribucnim centru j
    @variable(model, X, Bin) # BOOL 0 nebo 1, 1 - distr. centrum j se nachazi, 0 - jinak (nenachazi)
    @objective(model, Max, sum{P[s]*D[s][k]*W[k]*A[k]*f[s][j][k], s=1:@NUM_SCENARIOS, k=1:@NUM_ITEMS})
    @constraint()
    @constraint()
    @constraint()
    @constraint()
    @constraint()
    @constraint()
end
