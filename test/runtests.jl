using Random, GeometryBasics
using Query
using GeometryBasics: Polygon, MultiPolygon, Point

polys = [Polygon(rand(Point{2, Float32}, 20)) for i in 1:10]
pnames = [randstring(4) for i in 1:10]
numbers = rand(10)
bin = rand(Bool, 10)

x = MultiPolygon(polys, names = pnames, value = numbers, category = bin)

x = @from i in x begin
    @where i.value < 0.7
    @select i
    @collect
end
