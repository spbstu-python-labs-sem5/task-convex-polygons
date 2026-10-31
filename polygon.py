from __future__ import annotations

from math import fsum, hypot, isfinite, ulp
from typing import Iterable, Sequence

Point = tuple[float, float]
_EPS = 1e-9


def _cross(a: Point, b: Point, c: Point) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _coordinate_resolution(points: Iterable[Point]) -> float:
    return max((ulp(abs(value)) for point in points for value in point), default=0.0)


def _cross_tolerance(a: Point, b: Point, c: Point) -> float:
    scale = max(
        hypot(b[0] - a[0], b[1] - a[1]),
        hypot(c[0] - b[0], c[1] - b[1]),
        hypot(a[0] - c[0], a[1] - c[1]),
    )
    resolution = _coordinate_resolution((a, b, c))
    return _EPS * scale * scale + 4 * resolution * scale


def _area_tolerance(vertices: Sequence[Point]) -> float:
    origin = vertices[0]
    scale = max(hypot(x - origin[0], y - origin[1]) for x, y in vertices)
    resolution = _coordinate_resolution(vertices)
    return _EPS * scale * scale + 4 * resolution * scale


def _linear_tolerance(points: Sequence[Point]) -> float:
    if not points:
        return 0.0
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    scale = max(max(xs) - min(xs), max(ys) - min(ys))
    return _EPS * scale + 4 * _coordinate_resolution(points)


def _determinant_tolerance(
    first_length: float, second_length: float, points: Iterable[Point]
) -> float:
    resolution = _coordinate_resolution(points)
    return _EPS * first_length * second_length + 4 * resolution * (first_length + second_length)


class ConvexPolygon:
    """выпуклый многоугольник на плоскости."""

    def __init__(self, vertices: Iterable[Sequence[float]]) -> None:
        """создаёт многоугольник из вершин."""
        points: list[Point] = []
        for point in vertices:
            try:
                coordinates = tuple(point)
            except TypeError as error:
                raise ValueError("Каждая вершина должна содержать ровно две координаты") from error
            if len(coordinates) != 2:
                raise ValueError("Каждая вершина должна содержать ровно две координаты")
            try:
                converted = (float(coordinates[0]), float(coordinates[1]))
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError("Координаты вершин должны быть конечными числами") from error
            if not all(isfinite(value) for value in converted):
                raise ValueError("Координаты вершин должны быть конечными числами")
            points.append(converted)

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
        nonzero = [
            turn
            for i, turn in enumerate(turns)
            if abs(turn)
            > _cross_tolerance(points[i - 1], points[i], points[(i + 1) % len(points)])
        ]
        if not nonzero or (min(nonzero) < 0 < max(nonzero)):
            raise ValueError("Вершины должны задавать выпуклый многоугольник в порядке обхода границы")

        self._vertices = tuple(points)
        self._orientation = 1.0 if sum(nonzero) > 0 else -1.0
        for i, start in enumerate(points):
            end = points[(i + 1) % len(points)]
            if any(
                self._orientation * _cross(start, end, point)
                < -_cross_tolerance(start, end, point)
                for point in points
            ):
                raise ValueError("Вершины должны задавать выпуклый многоугольник в порядке обхода границы")
        if self.area <= _area_tolerance(self._vertices):
            raise ValueError("Площадь многоугольника должна быть положительной")

    @property
    def vertices(self) -> tuple[Point, ...]:
        """вершины по границе."""
        return self._vertices

    @property
    def area(self) -> float:
        """площадь многоугольника."""
        vertices = self._vertices
        origin_x, origin_y = vertices[0]
        return abs(
            fsum(
                (vertices[i][0] - origin_x) * (vertices[(i + 1) % len(vertices)][1] - origin_y)
                - (vertices[(i + 1) % len(vertices)][0] - origin_x) * (vertices[i][1] - origin_y)
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
            >= -_cross_tolerance(
                self._vertices[i], self._vertices[(i + 1) % len(self._vertices)], point
            )
            for i in range(len(self._vertices))
        )

    def triangulate(self) -> list[tuple[Point, Point, Point]]:
        """делит многоугольник на треугольники."""
        vertices = self._vertices
        return [
            (vertices[0], vertices[i], vertices[i + 1])
            for i in range(1, len(vertices) - 1)
            if abs(_cross(vertices[0], vertices[i], vertices[i + 1]))
            > _cross_tolerance(vertices[0], vertices[i], vertices[i + 1])
        ]

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
            previous_inside = (
                other._orientation * _cross(edge_start, edge_end, previous)
                >= -_cross_tolerance(edge_start, edge_end, previous)
            )
            for current in incoming:
                current_inside = (
                    other._orientation * _cross(edge_start, edge_end, current)
                    >= -_cross_tolerance(edge_start, edge_end, current)
                )
                if current_inside != previous_inside:
                    dx, dy = current[0] - previous[0], current[1] - previous[1]
                    ex, ey = edge_end[0] - edge_start[0], edge_end[1] - edge_start[1]
                    denominator = dx * ey - dy * ex
                    if abs(denominator) > _determinant_tolerance(
                        hypot(dx, dy), hypot(ex, ey), (previous, current, edge_start, edge_end)
                    ):
                        t = (
                            (edge_start[0] - previous[0]) * ey
                            - (edge_start[1] - previous[1]) * ex
                        ) / denominator
                        output.append((previous[0] + t * dx, previous[1] + t * dy))
                if current_inside:
                    output.append(current)
                previous, previous_inside = current, current_inside

        cleaned: list[Point] = []
        duplicate_tolerance = _linear_tolerance(output)
        for point in output:
            if (
                not cleaned
                or hypot(point[0] - cleaned[-1][0], point[1] - cleaned[-1][1])
                > duplicate_tolerance
            ):
                cleaned.append(point)
        if len(cleaned) > 1 and hypot(
            cleaned[0][0] - cleaned[-1][0], cleaned[0][1] - cleaned[-1][1]
        ) <= duplicate_tolerance:
            cleaned.pop()
        if len(cleaned) < 3:
            return None
        try:
            result = ConvexPolygon(cleaned)
        except ValueError:
            return None
        return result

    def __repr__(self) -> str:
        return f"ConvexPolygon({list(self._vertices)!r})"
