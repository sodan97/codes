/**
 * Carte stylisée de l'Afrique en matrice de points, calculée au build (aucune image externe).
 * Contour volontairement simplifié (longitude, latitude).
 */
type Point = [number, number];

const mainland: Point[] = [
  [-5.9, 35.8], [-1.0, 35.1], [3.0, 36.8], [10.2, 37.2], [11.1, 35.2], [10.0, 33.9], [11.5, 33.1],
  [15.2, 32.3], [18.5, 30.3], [20.1, 32.1], [23.0, 32.6], [25.1, 31.6], [29.9, 31.2], [32.3, 31.3],
  [32.6, 29.9], [34.0, 27.0], [35.5, 23.9], [37.2, 21.0], [38.6, 18.0], [39.5, 15.5], [41.7, 13.4],
  [43.2, 11.6], [44.0, 10.4], [51.2, 11.8], [51.0, 10.4], [49.5, 6.8], [48.0, 4.5], [45.3, 2.0],
  [42.5, -0.4], [39.7, -4.0], [39.3, -6.8], [40.5, -10.5], [40.7, -14.5], [37.0, -17.5], [34.8, -19.8],
  [35.5, -22.0], [35.4, -24.0], [32.6, -25.9], [31.0, -29.9], [28.0, -32.8], [25.6, -34.0], [20.0, -34.8],
  [18.4, -34.0], [17.9, -32.5], [16.5, -28.6], [15.1, -26.6], [14.5, -22.9], [11.8, -17.2], [13.4, -12.6],
  [13.2, -8.8], [12.2, -6.0], [11.8, -4.8], [8.7, -0.7], [9.4, 0.4], [9.8, 3.0], [9.7, 4.0], [8.5, 4.5],
  [6.0, 4.3], [4.5, 6.3], [3.4, 6.4], [-0.2, 5.6], [-2.0, 4.8], [-4.0, 5.3], [-7.5, 4.4], [-10.8, 6.3],
  [-13.2, 8.5], [-13.7, 9.5], [-15.0, 10.8], [-15.6, 11.9], [-16.8, 13.2], [-17.5, 14.7], [-16.5, 16.0],
  [-16.0, 18.1], [-17.05, 20.8], [-16.0, 23.7], [-14.5, 26.1], [-13.0, 27.6], [-11.0, 28.7], [-9.6, 30.4],
  [-9.8, 31.5], [-7.6, 33.6], [-6.8, 34.0],
];

const madagascar: Point[] = [
  [49.3, -12.0], [50.5, -15.5], [49.6, -17.0], [49.4, -18.0], [47.2, -25.0], [45.1, -25.5],
  [43.7, -23.5], [43.3, -22.0], [44.0, -17.0], [44.4, -16.2], [46.3, -15.7], [48.0, -13.5],
];

function inside([x, y]: Point, polygon: Point[]) {
  let hit = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
    const [xi, yi] = polygon[i];
    const [xj, yj] = polygon[j];
    if (yi > y !== yj > y && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) hit = !hit;
  }
  return hit;
}

export type MapDot = { x: number; y: number; region: "senegal" | "west" | "francophone" | "africa" };

const LON_MIN = -19;
const LON_MAX = 53;
const LAT_MAX = 38;
const LAT_MIN = -36;
export const MAP_WIDTH = 560;
export const MAP_HEIGHT = Math.round(((LAT_MAX - LAT_MIN) / (LON_MAX - LON_MIN)) * MAP_WIDTH);

export const project = (lon: number, lat: number) => ({
  x: Math.round(((lon - LON_MIN) / (LON_MAX - LON_MIN)) * MAP_WIDTH * 10) / 10,
  y: Math.round(((LAT_MAX - lat) / (LAT_MAX - LAT_MIN)) * MAP_HEIGHT * 10) / 10,
});

/** Zones très approximatives, utilisées uniquement pour la mise en valeur graphique. */
function regionOf(lon: number, lat: number): MapDot["region"] {
  if (lon >= -17.6 && lon <= -11.3 && lat >= 12.3 && lat <= 16.7) return "senegal";
  if (lon <= 16 && lat >= 4 && lat <= 25) return "west";
  const centralFrancophone = lon > 8 && lon <= 31 && lat >= -6 && lat < 24;
  const maghreb = lat >= 25 && lon <= 11.5;
  const indianOcean = lon >= 43 && lat <= -11;
  if (centralFrancophone || maghreb || indianOcean) return "francophone";
  return "africa";
}

export function africaDots(step = 1.6): MapDot[] {
  const dots: MapDot[] = [];
  for (let lat = LAT_MAX; lat >= LAT_MIN; lat -= step) {
    for (let lon = LON_MIN; lon <= LON_MAX; lon += step) {
      const point: Point = [lon, lat];
      if (inside(point, mainland) || inside(point, madagascar)) {
        dots.push({ ...project(lon, lat), region: regionOf(lon, lat) });
      }
    }
  }
  return dots;
}

export const DAKAR = project(-17.45, 14.7);
