# Generated by Django 4.2.5 on 2024-08-19 09:21

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("rolls", "0005_delete_analisi1_analisi_feature_alter_analisi_gene_and_more"),
    ]

    operations = [
        migrations.CreateModel(
            name="Gene",
            fields=[
                (
                    "id",
                    models.BigAutoField(
                        auto_created=True,
                        primary_key=True,
                        serialize=False,
                        verbose_name="ID",
                    ),
                ),
                ("name", models.CharField(max_length=100)),
            ],
        ),
        migrations.AddField(
            model_name="analisi",
            name="miRNA",
            field=models.CharField(blank=True, max_length=20, null=True),
        ),
        migrations.AlterField(
            model_name="analisi",
            name="feature",
            field=models.CharField(
                choices=[
                    (None, "Choice.."),
                    ("gender", "Gender"),
                    ("age_at_initial_pathologic_diagnosis", "Age"),
                    ("radiation_therapy", "Radiation therapy"),
                    ("patient_status", "Patient status"),
                    ("diabetes", "Diabetes"),
                    ("tobacco_smoking_history", "Tobacco smoking history"),
                    ("menopause_status", "Menopause status"),
                    ("alcohol_history_documented", "Alcohol history documented"),
                    ("pathologic_stage", "Pathologic stage"),
                ],
                default="Choice..",
                max_length=50,
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
                    ("LUSH", "LUSH"),
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
    ]
