from __future__ import annotations

from math import hypot
from typing import Iterable, Sequence

Point = tuple[float, float]
_EPS = 1e-9


def _cross(a: Point, b: Point, c: Point) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


class ConvexPolygon:
    """выпуклый многоугольник на плоскости."""

    def __init__(self, vertices: Iterable[Sequence[float]]) -> None:
        """создаёт многоугольник из вершин."""
        points = [(float(point[0]), float(point[1])) for point in vertices]
        if len(points) > 1 and points[0] == points[-1]:
            points.pop()
        if len(points) < 3:
            raise ValueError("У многоугольника должно быть не меньше трёх вершин")
        if len(set(points)) < 3:
            raise ValueError("У многоугольника должно быть не меньше трёх различных вершин")
        if any(points[i] == points[(i + 1) % len(points)] for i in range(len(points))):
            raise ValueError("Соседние вершины должны различаться")

        turns = [
            _cross(points[i - 1], points[i], points[(i + 1) % len(points)])
            for i in range(len(points))
        ]
        nonzero = [turn for turn in turns if abs(turn) > _EPS]
        if not nonzero or (min(nonzero) < 0 < max(nonzero)):
            raise ValueError("Вершины должны задавать выпуклый многоугольник в порядке обхода границы")

        self._vertices = tuple(points)
        self._orientation = 1.0 if sum(nonzero) > 0 else -1.0
        for i, start in enumerate(points):
            end = points[(i + 1) % len(points)]
            if any(self._orientation * _cross(start, end, point) < -_EPS for point in points):
                raise ValueError("Вершины должны задавать выпуклый многоугольник в порядке обхода границы")
        if self.area <= _EPS:
            raise ValueError("Площадь многоугольника должна быть положительной")

    @property
    def vertices(self) -> tuple[Point, ...]:
        """вершины по границе."""
        return self._vertices

    @property
    def area(self) -> float:
        """площадь многоугольника."""
        vertices = self._vertices
        return abs(
            sum(
                vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
                - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
                for i in range(len(vertices))
            ) / 2
        )

    @property
    def perimeter(self) -> float:
        """периметр многоугольника."""
        vertices = self._vertices
        return sum(
            hypot(
                vertices[i][0] - vertices[(i + 1) % len(vertices)][0],
                vertices[i][1] - vertices[(i + 1) % len(vertices)][1],
            )
            for i in range(len(vertices))
        )

    def contains(self, point: Sequence[float]) -> bool:
        """проверяет, внутри ли точка или на границе."""
        point = (float(point[0]), float(point[1]))
        return all(
            self._orientation
            * _cross(self._vertices[i], self._vertices[(i + 1) % len(self._vertices)], point)
            >= -_EPS
            for i in range(len(self._vertices))
        )

    def triangulate(self) -> list[tuple[Point, Point, Point]]:
        """делит многоугольник на треугольники."""
        vertices = self._vertices
        return [(vertices[0], vertices[i], vertices[i + 1]) for i in range(1, len(vertices) - 1)]

    def intersection(self, other: ConvexPolygon) -> ConvexPolygon | None:
        """находит пересечение с другим многоугольником."""
        output = list(self._vertices)
        clip = other._vertices
        for i, edge_start in enumerate(clip):
            edge_end = clip[(i + 1) % len(clip)]
            incoming = output
            output = []
            if not incoming:
                return None
            previous = incoming[-1]
            previous_inside = other._orientation * _cross(edge_start, edge_end, previous) >= -_EPS
            for current in incoming:
                current_inside = other._orientation * _cross(edge_start, edge_end, current) >= -_EPS
                if current_inside != previous_inside:
                    dx, dy = current[0] - previous[0], current[1] - previous[1]
                    ex, ey = edge_end[0] - edge_start[0], edge_end[1] - edge_start[1]
                    denominator = dx * ey - dy * ex
                    if abs(denominator) > _EPS:
                        t = (
                            (edge_start[0] - previous[0]) * ey
                            - (edge_start[1] - previous[1]) * ex
                        ) / denominator
                        output.append((previous[0] + t * dx, previous[1] + t * dy))
                if current_inside:
                    output.append(current)
                previous, previous_inside = current, current_inside

        cleaned: list[Point] = []
        for point in output:
            if not cleaned or hypot(point[0] - cleaned[-1][0], point[1] - cleaned[-1][1]) > _EPS:
                cleaned.append(point)
        if len(cleaned) > 1 and hypot(
            cleaned[0][0] - cleaned[-1][0], cleaned[0][1] - cleaned[-1][1]
        ) <= _EPS:
            cleaned.pop()
        if len(cleaned) < 3:
            return None
        try:
            result = ConvexPolygon(cleaned)
        except ValueError:
            return None
        return result if result.area > _EPS else None

    def __repr__(self) -> str:
        return f"ConvexPolygon({list(self._vertices)!r})"
