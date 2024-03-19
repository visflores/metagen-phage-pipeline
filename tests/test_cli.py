import shutil
import os
import subprocess
import unittest
from pathlib import Path

class MetagenTestCli(unittest.TestCase):
    """
        Classe para teste simples de funcionamento da pipeline.

        Esse teste pode demorar alguns minutos.
    """

    BASE = Path(__file__).parent

    def tearDown(self):
        """
            Criação de diretório temporário para resultados e incluindo ele
            para remoção após testes.
        """
        results = MetagenTestCli.BASE / "results"
        input = MetagenTestCli.BASE / "input.json"
        output = MetagenTestCli.BASE / "output.json"
        r1 = MetagenTestCli.BASE / "assets" / "reads" / "multi_1.fastq"
        r2 = MetagenTestCli.BASE / "assets" / "reads" / "multi_2.fastq"

        self.addCleanup(shutil.rmtree, results)
        self.addCleanup(os.remove, input)
        self.addCleanup(os.remove, output)
        self.addCleanup(os.remove, r1)
        self.addCleanup(os.remove, r2)

    def _runMetagen(self, multifasta, reads):
        """
            Método para executar a pipeline.
        """
        cmd = [
            "metagen-phage",
            "-file",
            str(multifasta),
            "-reads",
            str(reads),
            "-t",
            "4"
        ]

        run = subprocess.run(cmd, capture_output=True)
        return run.returncode

    def test_metagen(self):
        """
            Executa teste.
        """
        multifasta = Path(__file__).parent / "assets" / "multifasta.fasta"
        reads = Path(__file__).parent / "assets" / "reads"

        result = self._runMetagen(multifasta, reads)

        self.assertEqual(result, 0)



