import { defineConfig } from 'vite';

// Relative asset paths: the same build is served by Tauri and by the
// development bridge (examples/editor_server.rs).
export default defineConfig({
  base: './',
  build: { outDir: 'dist', emptyOutDir: true, chunkSizeWarningLimit: 2000 },
  server: { proxy: { '/api': 'http://127.0.0.1:8787' } },
});
