"""
    O objetivo desse script é o de criar uma classe para auxliar a interface de
    linha de comando com a execução da pipeline de isolamento de bacteriófagos.

    Caso queira adicionar novas opções à pipeline, lembre-se de que os métodos
    para lidar com essas novas opções devem ser depositados nessa classe.
    Portanto, qualquer extensão da pipeline que exija novos inputs e que
    serão incorporados a CLI devem refletirem também nessa classe.
"""
import subprocess
import json

from pathlib import Path


class CliAux:
    """
        Classe que auxilia a CLI para montagem do JSON de input e execução da
        pipeline de recuperação de bacteriófagos.

        Author
        ------
            Vinicius S. Flores <vinicius.nuvem@gmail.com>

    """
    BASE_PATH = Path(__file__).parents[1]
    CONFIG_PATH = BASE_PATH / "pipes" / "config" / "config.cfg"
    FLOW = BASE_PATH / "pipes" / "workflows" / "main_flow.wdl"

    def __init__(self,
                 scaffold_file,
                 reads_dir,
                 threads=4,
                 find_prophages=False
                ):
        """
            Construtor da classe CliAux.

            Parameters
            ----------
            scaffold_file : Path
                Caminho do arquivo multifasta contendo os contigs que foram
                montados a partir dos reads passados para o parâmetro
                'reads_dir'. Esse arquivo não deve estar comprimido (gziped).

            reads_dir : Path
                Caminho do diretório contendo os reads R1 e R2 usados para a
                construção dos contigs presentes no arquivo multifasta passado
                para o parâmetro 'scaffold_file'.

            threads : int
                Quantidade de threads que será usada pelas ferramentas
                presentes na pipeline.

            find_prophages : bool
                Ativa ou desativa a busca por profagos nos MAGs montados.
        """
        self._scaffold_file = scaffold_file
        self._reads_dir = reads_dir
        self._threads = threads
        self._find_prophages = find_prophages

    def exec_pipe(self):
        self.__create_json()

        cmd = [
            "miniwdl",
            "run",
            "-i",
            self._input_json.absolute().as_posix(),
            "--cfg",
            CliAux.CONFIG_PATH.absolute().as_posix(),
            "-o",
            "output.json",
            "-d",
            "./results/",
            CliAux.FLOW.absolute().as_posix()
        ]

        runned = subprocess.run(cmd)

        if runned.returncode != 0:
            raise CalledProcessError("Erro ao executar pipeline!")

    def __create_json(self):
        """
            Cria o JSON de input da pipeline de predição de bacteriófagos
        """
        self.__group_reads()

        dicio = {}

        dicio["MainFlow.multiFasta"] = self.scaffold_file.absolute().as_posix()
        dicio["MainFlow.readsToMap"] = self.grouped_reads
        dicio["MainFlow.threads"] = self.threads
        dicio["MainFlow.findProvirus"] = self.find_prophages

        with open("input.json", "w") as input:
            json.dump(dicio, input, indent=4)

        self._input_json = Path("input.json")


    def __group_reads(self):
        """
            Método para agrupamento dos reads R1 e R2 e formação da estrutura
            usada pela pipeline.

        Return
        ------
        self.__grouped_reads : List[Dict]
            Lista contendo dicionários com os respectivos reads R1 e R2
            agrupados.
        """
        r1, r2 = (
            list(self.reads_dir.glob("*_1.fastq")),
            list(self.reads_dir.glob("*_2.fastq"))
        )

        grouped_reads = []

        for one, two in zip(r1, r2):
            dicio = {
                "left": one.absolute().as_posix(),
                "right": two.absolute().as_posix()
            }

            grouped_reads.append(dicio)

        self.__grouped_reads = grouped_reads


    @property
    def grouped_reads(self):
        return self.__grouped_reads

    @property
    def scaffold_file(self):
        return self._scaffold_file

    @property
    def reads_dir(self):
        return self._reads_dir

    @property
    def threads(self):
        return self._threads

    @property
    def find_prophages(self):
        return self._find_prophages

