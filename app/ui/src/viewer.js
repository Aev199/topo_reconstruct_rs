// 3D view of the reconstructed geometry: surfaces, contour and junction
// edges, bars, audit findings, selection and picking.
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

const CLASS_COLOR = { failure: 0xd92d20, plaxis: 0xe07b00, review: 0x7a8394 };

function stiffnessColor(k) {
  // Stable, well separated hues per stiffness id.
  const hue = ((k * 0.61803398875) % 1 + 1) % 1;
  return new THREE.Color().setHSL(hue, 0.35, 0.72);
}

export class Viewer {
  constructor(container, onPick) {
    this.container = container;
    this.onPick = onPick;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(window.devicePixelRatio || 1);
    this.renderer.setClearColor(0xeef0f4);
    container.appendChild(this.renderer.domElement);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1e6);
    this.camera.up.set(0, 0, 1);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = false;
    this.controls.addEventListener('change', () => this.render());
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x8890a0, 1.6));
    const sun = new THREE.DirectionalLight(0xffffff, 1.2);
    sun.position.set(1, 0.6, 2);
    this.scene.add(sun);
    this.model = new THREE.Group();
    this.marks = new THREE.Group();
    this.selected = new THREE.Group();
    this.scene.add(this.model, this.marks, this.selected);
    this.raycaster = new THREE.Raycaster();
    this.offset = new THREE.Vector3();
    this.data = null;
    this.size = 1;
    let down = null;
    const canvas = this.renderer.domElement;
    canvas.addEventListener('pointerdown', (e) => { down = [e.clientX, e.clientY]; });
    canvas.addEventListener('pointerup', (e) => {
      if (down && Math.hypot(e.clientX - down[0], e.clientY - down[1]) < 5) this.onPick(e);
      down = null;
    });
    new ResizeObserver(() => this.resize()).observe(container);
    this.resize();
  }

  resize() {
    const { clientWidth: w, clientHeight: h } = this.container;
    if (!w || !h) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.render();
  }

  render() {
    this.renderer.render(this.scene, this.camera);
  }

  /// Local coordinates: the model is re-centred to keep float32 precision.
  local(p) {
    return new THREE.Vector3(p[0] - this.offset.x, p[1] - this.offset.y, p[2] - this.offset.z);
  }

  world(v) {
    return [v.x + this.offset.x, v.y + this.offset.y, v.z + this.offset.z];
  }

  setScene(data) {
    const first = !this.data;
    this.data = data;
    const [lo, hi] = data.bounds;
    if (first) {
      this.offset.set((lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2);
    }
    this.size = Math.max(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-3);
    for (const child of [...this.model.children]) {
      child.geometry.dispose();
      child.material.dispose();
      this.model.remove(child);
    }
    const local = data.vertices.map((p) => [p[0] - this.offset.x, p[1] - this.offset.y, p[2] - this.offset.z]);
    this.localVertices = local;

    // Surfaces: one geometry, a triangle-to-surface map for picking.
    let count = 0;
    for (const s of data.surfaces) count += s.triangles.length;
    const positions = new Float32Array(count * 3);
    const colors = new Float32Array(count * 3);
    this.triangleSurface = new Int32Array(count / 3);
    this.triangleVertices = new Int32Array(count);
    let k = 0;
    data.surfaces.forEach((s, index) => {
      const c = stiffnessColor(s.stiffness);
      for (let i = 0; i < s.triangles.length; i++) {
        const p = local[s.triangles[i]];
        positions.set(p, k * 3);
        colors.set([c.r, c.g, c.b], k * 3);
        this.triangleVertices[k] = s.triangles[i];
        if (i % 3 === 0) this.triangleSurface[k / 3] = index;
        k++;
      }
    });
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.computeVertexNormals();
    this.surfaceMesh = new THREE.Mesh(geometry, new THREE.MeshLambertMaterial({
      vertexColors: true, side: THREE.DoubleSide, polygonOffset: true, polygonOffsetFactor: 1, polygonOffsetUnits: 1,
    }));
    this.model.add(this.surfaceMesh);

    // Edges: contour (dark) and embedded junction lines (blue).
    const edgePositions = new Float32Array(data.edges.length * 6);
    const edgeColors = new Float32Array(data.edges.length * 6);
    data.edges.forEach(([, a, b, embedded], i) => {
      edgePositions.set(local[a], i * 6);
      edgePositions.set(local[b], i * 6 + 3);
      const c = embedded ? [0.15, 0.35, 0.85] : [0.2, 0.22, 0.27];
      edgeColors.set(c, i * 6);
      edgeColors.set(c, i * 6 + 3);
    });
    const edgeGeometry = new THREE.BufferGeometry();
    edgeGeometry.setAttribute('position', new THREE.BufferAttribute(edgePositions, 3));
    edgeGeometry.setAttribute('color', new THREE.BufferAttribute(edgeColors, 3));
    this.edgeLines = new THREE.LineSegments(edgeGeometry, new THREE.LineBasicMaterial({ vertexColors: true }));
    this.model.add(this.edgeLines);

    // Bars.
    const barPositions = new Float32Array(data.bars.length * 6);
    data.bars.forEach(([, a, b], i) => {
      barPositions.set(local[a], i * 6);
      barPositions.set(local[b], i * 6 + 3);
    });
    const barGeometry = new THREE.BufferGeometry();
    barGeometry.setAttribute('position', new THREE.BufferAttribute(barPositions, 3));
    this.barLines = new THREE.LineSegments(barGeometry, new THREE.LineBasicMaterial({ color: 0xc2410c }));
    this.model.add(this.barLines);

    // Vertices in use, on a grid for nearest-vertex picking.
    const used = new Set();
    for (const [, a, b] of data.edges) { used.add(a); used.add(b); }
    for (const [, a, b] of data.bars) { used.add(a); used.add(b); }
    this.cell = this.size / 200;
    this.grid = new Map();
    for (const v of used) {
      const key = this.key(local[v]);
      if (!this.grid.has(key)) this.grid.set(key, []);
      this.grid.get(key).push(v);
    }
    if (first) this.fit();
    this.render();
  }

  key(p) {
    return p.map((x) => Math.floor(x / this.cell)).join(',');
  }

  nearestVertex(point, radius) {
    const c = [point.x, point.y, point.z].map((x) => Math.floor(x / this.cell));
    const r = Math.max(1, Math.ceil(radius / this.cell));
    let best = null;
    let bestD = radius;
    for (let x = -r; x <= r; x++) for (let y = -r; y <= r; y++) for (let z = -r; z <= r; z++) {
      for (const v of this.grid.get(`${c[0] + x},${c[1] + y},${c[2] + z}`) || []) {
        const p = this.localVertices[v];
        const d = Math.hypot(p[0] - point.x, p[1] - point.y, p[2] - point.z);
        if (d <= bestD) { bestD = d; best = v; }
      }
    }
    return best;
  }

  fit() {
    const d = this.size * 1.4;
    this.controls.target.set(0, 0, 0);
    this.camera.position.set(d * 0.8, -d, d * 0.7);
    this.camera.near = this.size / 1e4;
    this.camera.far = this.size * 100;
    this.camera.updateProjectionMatrix();
    this.controls.update();
    this.render();
  }

  flyTo(point, extent) {
    const p = this.local(point);
    const d = Math.max(extent || 0, this.size / 50, 0.5) * 3;
    const dir = this.camera.position.clone().sub(this.controls.target).normalize();
    this.controls.target.copy(p);
    this.camera.position.copy(p.clone().add(dir.multiplyScalar(d)));
    this.controls.update();
    this.render();
  }

  /// Audit findings as markers: points (pixel size) and segments.
  setFindings(findings) {
    for (const child of [...this.marks.children]) {
      child.geometry.dispose();
      child.material.dispose();
      this.marks.remove(child);
    }
    for (const cls of ['review', 'plaxis', 'failure']) {
      const list = findings.filter((f) => f.class === cls);
      if (!list.length) continue;
      const points = list.map((f) => this.local(f.points[0]));
      const pg = new THREE.BufferGeometry().setFromPoints(points);
      this.marks.add(new THREE.Points(pg, new THREE.PointsMaterial({
        color: CLASS_COLOR[cls], size: cls === 'failure' ? 11 : 8, sizeAttenuation: false, depthTest: false,
      })));
      const segments = list.filter((f) => f.points.length > 1)
        .flatMap((f) => [this.local(f.points[0]), this.local(f.points[1])]);
      if (segments.length) {
        const sg = new THREE.BufferGeometry().setFromPoints(segments);
        this.marks.add(new THREE.LineSegments(sg, new THREE.LineBasicMaterial({ color: CLASS_COLOR[cls], depthTest: false })));
      }
    }
    this.marks.renderOrder = 2;
    this.render();
  }

  /// Highlight a selection: surfaces, vertices, edges, bars, a finding.
  highlight(sel) {
    for (const child of [...this.selected.children]) {
      child.geometry.dispose();
      child.material.dispose();
      this.selected.remove(child);
    }
    if (!this.data || !sel) return this.render();
    const yellow = 0xfacc15;
    for (const s of sel.surfaces || []) {
      const tri = this.data.surfaces[s]?.triangles || [];
      const g = new THREE.BufferGeometry().setFromPoints(tri.map((v) => new THREE.Vector3(...this.localVertices[v])));
      this.selected.add(new THREE.Mesh(g, new THREE.MeshBasicMaterial({
        color: yellow, transparent: true, opacity: 0.55, side: THREE.DoubleSide, depthTest: false,
      })));
    }
    const pts = (sel.vertices || []).map((v) => new THREE.Vector3(...this.localVertices[v]));
    if (sel.point) pts.push(this.local(sel.point));
    if (pts.length) {
      this.selected.add(new THREE.Points(new THREE.BufferGeometry().setFromPoints(pts),
        new THREE.PointsMaterial({ color: 0x111827, size: 12, sizeAttenuation: false, depthTest: false })));
    }
    const segs = [];
    for (const e of sel.edges || []) {
      const edge = this.data.edges.find((x) => x[0] === e);
      if (edge) segs.push(new THREE.Vector3(...this.localVertices[edge[1]]), new THREE.Vector3(...this.localVertices[edge[2]]));
    }
    for (const b of sel.bars || []) {
      for (const [axis, u, w] of this.data.bars) {
        if (axis === b) segs.push(new THREE.Vector3(...this.localVertices[u]), new THREE.Vector3(...this.localVertices[w]));
      }
    }
    if (sel.segment) segs.push(this.local(sel.segment[0]), this.local(sel.segment[1]));
    if (segs.length) {
      this.selected.add(new THREE.LineSegments(new THREE.BufferGeometry().setFromPoints(segs),
        new THREE.LineBasicMaterial({ color: 0x111827, depthTest: false })));
    }
    this.selected.renderOrder = 3;
    this.render();
  }

  /// What is under the pointer in a pick mode.
  pick(event, mode) {
    if (!this.data) return null;
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(((event.clientX - rect.left) / rect.width) * 2 - 1,
      -((event.clientY - rect.top) / rect.height) * 2 + 1);
    this.raycaster.setFromCamera(ndc, this.camera);
    const distance = this.camera.position.distanceTo(this.controls.target);
    const tolerance = distance * 0.008;
    this.raycaster.params.Line.threshold = tolerance;
    if (mode === 'edge' || mode === 'bar') {
      const target = mode === 'edge' ? this.edgeLines : this.barLines;
      const hit = this.raycaster.intersectObject(target)[0];
      if (!hit) return null;
      const segment = Math.floor(hit.index / 2);
      const point = this.world(hit.point);
      return mode === 'edge'
        ? { edge: this.data.edges[segment][0], point }
        : { bar: this.data.bars[segment][0], point };
    }
    const hits = this.raycaster.intersectObjects([this.surfaceMesh, this.edgeLines, this.barLines]);
    if (!hits.length) return null;
    const hit = hits[0];
    if (mode === 'vertex') {
      const v = this.nearestVertex(hit.point, tolerance * 2);
      return v === null ? null : { vertex: v, point: this.data.vertices[v] };
    }
    const surfaceHit = hits.find((h) => h.object === this.surfaceMesh);
    if (!surfaceHit) return null;
    return { surface: this.triangleSurface[surfaceHit.faceIndex], point: this.world(surfaceHit.point) };
  }
}
