import subprocess
import sys
import unittest

HEAVY_DEPS = ["numpy", "pandas", "matplotlib", "scipy", "numba", "plotly"]


class TestStartupImports(unittest.TestCase):
    def test_no_heavy_deps_at_startup(self):
        """Importing seismicrna.main must not load any heavy dependencies."""
        script = (
            "import sys;"
            "import seismicrna.main;"
            "heavy = [m for m in sys.modules"
            f" if any(m == d or m.startswith(d + '.') for d in {HEAVY_DEPS!r})];"
            "print('\\n'.join(heavy))"
        )
        result = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        loaded = [m for m in result.stdout.splitlines() if m]
        self.assertEqual(
            loaded, [], msg=f"Heavy dependencies loaded at startup: {loaded}"
        )
