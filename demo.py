from polygon import ConvexPolygon

square = ConvexPolygon([(0, 0), (4, 0), (4, 4), (0, 4)])
triangle = ConvexPolygon([(2, 1), (5, 1), (2, 4)])
intersection = square.intersection(triangle)

print("Площадь квадрата:", square.area)
print("Периметр квадрата:", square.perimeter)
print("Точка (1, 1) внутри квадрата:", "да" if square.contains((1, 1)) else "нет")
print("Точка (5, 5) внутри квадрата:", "да" if square.contains((5, 5)) else "нет")
print("Площадь пересечения:", intersection.area if intersection else 0)
print("Треугольники после триангуляции:")
for vertices in square.triangulate():
    print(vertices)

try:
    ConvexPolygon([(0, 0), (2, 0), (1, 1), (2, 2), (0, 2)])
except ValueError as error:
    print("Ошибка при проверке выпуклости:", error)
