"""Console script for soil_heat."""
import soil_heat

import typer
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def main():
    """Console script for soil_heat."""
    console.print("Replace this message by putting your code into "
               "soil_heat.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    


if __name__ == "__main__":
    app()
