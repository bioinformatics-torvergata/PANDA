# Generated by Django 5.1 on 2024-09-27 15:27

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('rolls', '0010_analisi_mutation'),
    ]

    operations = [
        migrations.CreateModel(
            name='Protein',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('protein', models.CharField(max_length=100)),
            ],
        ),
    ]
