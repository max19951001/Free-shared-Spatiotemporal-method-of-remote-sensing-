FUNCTION Cost_Fun, g
  common V_PUB1
  common V_PUB2
  L=total((ind_v ## g-dep_v)^2)
  RETURN, L
END