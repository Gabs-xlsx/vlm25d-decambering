# VLM 2,5D com *decambering*

Estimativa de cargas e momento de dobradiça em superfícies de comando de aeronaves em escala reduzida, acoplando o Método de Lençóis de Vórtices a polares bidimensionais viscosas (Mukherjee e Gopalarathnam, 2006). Código do TCC de Gabriel Luiz Santos Resende Gomes (Engenharia Aeroespacial, UFMG, 2026).

## Requisitos

- MATLAB R2020a ou mais recente.
- [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/): coloque o `xfoil.exe` nesta pasta.

## Como rodar

Abra esta pasta no MATLAB e execute:

```matlab
exemplo
```

O exemplo calcula a asa retangular de alongamento 6 com perfil S1223 (α = 5°, aileron em +10°) e imprime C_L e C_h das soluções potencial e acoplada. Para outro caso, altere a geometria, o perfil e a condição de voo no início do `exemplo.m`.

## Conteúdo

- `exemplo.m`: caso de uso completo.
- `src/`: rotinas da ferramenta.
- `Perfis/`: coordenadas dos perfis (base de dados da UIUC).
