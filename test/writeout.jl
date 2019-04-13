#Input: Mcormick object display from REPL
#   MC{3}(0.0, 0.0, [0, 0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], true)
#Output: How you would format the initialization of that object
#   MC{3}(0.0, 0.0, IntervalType(0, 0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), true)


function writeout(MC::String)
#Pieces
N = MC[4]
sv = "SVector{"*"$(N)"*",Float64}("
intv = "IntervalType("
MC = replace(MC, '['=>intv, count=1)
MC = replace(MC, '['=>sv, count=2)
MC = replace(MC, ']'=>')')

println(MC)
end

function writeout(MC2::MC) #Requires EAGO or some handler of ::MC to load
#Pieces
MC = "$(MC2)"
N = MC[4]
sv = "SVector{"*"$(N)"*",Float64}("
intv = "IntervalType("
MC = replace(MC, '['=>intv, count=1)
MC = replace(MC, '['=>sv, count=2)
MC = replace(MC, ']'=>')')

println(MC)
end
