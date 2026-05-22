import inspect
import unittest as ut

from click.testing import CliRunner

from seismicrna.core.logs import Level, get_config, set_config
from seismicrna.sim import (abstract, clusts, ends, fastq, fold,
                            muts, params, ref, idmut, total)


class TestSimCLIParams(ut.TestCase):
    """Check every CLI param name appears in run()'s signature.

    inspect.signature() follows the __wrapped__ chain (all decorators in
    run_func use @wraps) and stops at the __signature__ set by
    param_defaults, which preserves the original parameter names.
    """

    # Params consumed by decorator infrastructure before reaching run().
    _consumed_by_decorator = frozenset({"tmp_pfx"})

    def _assert_params_match_run(self, module):
        sig = inspect.signature(module.run)
        sig_param_names = set(sig.parameters)
        for cli_param in module.params:
            if cli_param.name in self._consumed_by_decorator:
                continue
            self.assertIn(
                cli_param.name, sig_param_names,
                f"[{module.COMMAND}] CLI param '{cli_param.name}' is in "
                f"params but missing from run()"
            )

    def test_abstract(self):  self._assert_params_match_run(abstract)
    def test_clusts(self):    self._assert_params_match_run(clusts)
    def test_ends(self):      self._assert_params_match_run(ends)
    def test_fastq(self):     self._assert_params_match_run(fastq)
    def test_fold(self):      self._assert_params_match_run(fold)
    def test_muts(self):      self._assert_params_match_run(muts)
    def test_params(self):    self._assert_params_match_run(params)
    def test_ref(self):       self._assert_params_match_run(ref)
    def test_idmut(self):    self._assert_params_match_run(idmut)
    def test_total(self):     self._assert_params_match_run(total)


class TestSimCLIInvocation(ut.TestCase):
    """Smoke-test each subcommand via CliRunner.

    Two tiers:
    - --help: validates Click's param parsing and help-text generation
      for every command (including fold, which needs a positional arg).
    - empty invocation ([]): calls run() through the full decorator
      stack with no files; commands return empty results with exit_code 0.
      fold is excluded because its positional arg is required and Click
      itself rejects an empty invocation with exit_code 2.
    """

    def setUp(self):
        self._config = get_config()
        set_config(verbosity=Level.ERROR, log_file_path=None,
                   exit_on_error=True)

    def tearDown(self):
        set_config(**self._config._asdict())

    def _invoke_help(self, module):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(module.cli, ["--help"])
        self.assertEqual(result.exit_code, 0, msg=result.output)

    def _invoke_empty(self, module):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(module.cli, [])
        self.assertEqual(result.exit_code, 0, msg=result.output)

    # --help tests (all commands)
    def test_abstract_help(self):  self._invoke_help(abstract)
    def test_clusts_help(self):    self._invoke_help(clusts)
    def test_ends_help(self):      self._invoke_help(ends)
    def test_fastq_help(self):     self._invoke_help(fastq)
    def test_fold_help(self):      self._invoke_help(fold)
    def test_muts_help(self):      self._invoke_help(muts)
    def test_params_help(self):    self._invoke_help(params)
    def test_ref_help(self):       self._invoke_help(ref)
    def test_idmut_help(self):    self._invoke_help(idmut)
    def test_total_help(self):     self._invoke_help(total)

    # Empty-invocation tests (all commands except fold)
    def test_abstract_empty(self): self._invoke_empty(abstract)
    def test_clusts_empty(self):   self._invoke_empty(clusts)
    def test_ends_empty(self):     self._invoke_empty(ends)
    def test_fastq_empty(self):    self._invoke_empty(fastq)
    def test_muts_empty(self):     self._invoke_empty(muts)
    def test_params_empty(self):   self._invoke_empty(params)
    def test_ref_empty(self):      self._invoke_empty(ref)
    def test_idmut_empty(self):   self._invoke_empty(idmut)
    def test_total_empty(self):    self._invoke_empty(total)


if __name__ == "__main__":
    ut.main()
