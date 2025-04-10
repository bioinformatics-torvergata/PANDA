# Generated by Django 4.2.5 on 2024-11-11 14:24

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("rolls", "0012_gene_symbol"),
    ]

    operations = [
        migrations.AddField(
            model_name="analisi_mutation",
            name="feature_selected",
            field=models.CharField(default="Choice..", max_length=50),
        ),
        migrations.AddField(
            model_name="analisi_mutation",
            name="number",
            field=models.IntegerField(
                choices=[(10, "10"), (15, "15"), (20, "20"), (25, "25")], default=10
            ),
        ),
        migrations.AlterField(
            model_name="analisi",
            name="tumor",
            field=models.CharField(
                choices=[
                    (None, "Choice.."),
                    ("ACC", "ACC"),
                    ("BLCA", "BLCA"),
                    ("BRCA", "BRCA"),
                    ("CESC", "CESC"),
                    ("CHOL", "CHOL"),
                    ("COAD", "COAD"),
                    ("DLBC", "DLBC"),
                    ("ESCA", "ESCA"),
                    ("GBM", "GBM"),
                    ("HNSC", "HNSC"),
                    ("KICH", "KICH"),
                    ("KIRC", "KIRC"),
                    ("KIRP", "KIRP"),
                    ("LGG", "LGG"),
                    ("LIHC", "LIHC"),
                    ("LUAD", "LUAD"),
                    ("LUSC", "LUSC"),
                    ("MESO", "MESO"),
                    ("OV", "OV"),
                    ("PAAD", "PAAD"),
                    ("PCPG", "PCPG"),
                    ("PRAD", "PRAD"),
                    ("READ", "READ"),
                    ("SARC", "SARC"),
                    ("SKCM", "SKCM"),
                    ("STAD", "STAD"),
                    ("TGCT", "TGCT"),
                    ("THCA", "THCA"),
                    ("THYM", "THYM"),
                    ("UCEC", "UCEC"),
                    ("UCS", "UCS"),
                    ("UVM", "UVM"),
                ],
                default="Choice..",
                max_length=10,
            ),
        ),
        migrations.AlterField(
            model_name="analisi_mutation",
            name="tumor",
            field=models.CharField(
                choices=[
                    (None, "Choice.."),
                    ("ACC", "ACC"),
                    ("BLCA", "BLCA"),
                    ("BRCA", "BRCA"),
                    ("CESC", "CESC"),
                    ("CHOL", "CHOL"),
                    ("COAD", "COAD"),
                    ("DLBC", "DLBC"),
                    ("ESCA", "ESCA"),
                    ("GBM", "GBM"),
                    ("HNSC", "HNSC"),
                    ("KICH", "KICH"),
                    ("KIRC", "KIRC"),
                    ("KIRP", "KIRP"),
                    ("LGG", "LGG"),
                    ("LIHC", "LIHC"),
                    ("LUAD", "LUAD"),
                    ("LUSC", "LUSC"),
                    ("OV", "OV"),
                    ("PAAD", "PAAD"),
                    ("PCPG", "PCPG"),
                    ("PRAD", "PRAD"),
                    ("READ", "READ"),
                    ("SARC", "SARC"),
                    ("SKCM", "SKCM"),
                    ("STAD", "STAD"),
                    ("TGCT", "TGCT"),
                    ("THCA", "THCA"),
                    ("THYM", "THYM"),
                    ("UCEC", "UCEC"),
                    ("UCS", "UCS"),
                    ("UVM", "UVM"),
                ],
                default="Choice..",
                max_length=10,
            ),
        ),
    ]
