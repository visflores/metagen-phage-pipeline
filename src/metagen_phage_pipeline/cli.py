import argparse

from pathlib import Path
from cli_aux.cli_aux import CliAux

def main():
    parser = argparse.ArgumentParser(
        prog="metage-phage-pipe",
        description="""Pipeline para a recuperação de MAGs de bacteriófagos a
        partir de metagenomas."""
    )

    parser.add_argument("-file",
                        "--scaffold-file",
                        required=True,
                        type=Path,
                        help="""Arquivo multifasta com todos os contigs
                        montados a partir dos reads passados a flag
                        '--reads-dir'."""
                       )

    parser.add_argument("-reads",
                        "--reads-dir",
                        required=True,
                        type=Path,
                        help=""""Diretório contendo os reads usados para a
                        montagem dos contigos do arquivo multifasta passado
                        para a flag '--scaffold-file'."""
                       )

    parser.add_argument("-t",
                        "--threads",
                        default=4,
                        required=False,
                        type=int,
                        help="""Quantidade de threads que será usada para a
                        execução da pipeline."""
                       )

    parser.add_argument("-pro",
                        "--prophage-find",
                        action="store_true",
                        default=False,
                        required=False,
                        help="""Ativa a busca por profagos nos MAGs montados.
                        Ativar essa opção torna o processo mais demorado."""
                       )

    args = parser.parse_args()

    obj = CliAux(args.scaffold_file,
                 args.reads_dir,
                 args.threads,
                 args.prophage_find)

    obj.exec_pipe()


if __name__ == "__main__":
    main()

