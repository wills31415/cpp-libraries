# Anciennes bibliothèques, en attente de refonte

Ces bibliothèques datent d'avant la refonte du dépôt. Elles sont **figées** : ni
construites, ni testées, ni installées. Elles restent lisibles comme références, et chacune
en sortira refondue selon la charte (design/interfaces.md, décision 36) quand son domaine
d'accueil sera conçu.

| Bibliothèque | Destination |
|---|---|
| `array2d` | remplacée par `wsl::mdarray` ; ses vues de lignes et de colonnes servent de référence au lot 2, puis elle disparaît |
| `graded`, `graded_io` | un scalaire étendu (jet en un infinitésimal positif, groupe ordonné non archimédien), avec la question des scalaires étendus ; renommée à ce moment-là ; `graded_io` rejoint la notion de son type |
| `elo` | `statistics`, catalogue de modèles de classement (Elo, Bradley-Terry) |
| `progress` | `execution`, suivi d'un calcul long |

Pour les lire avec leurs tests, telles qu'à la version `v0.1.0` : `git checkout v0.1.0`.
