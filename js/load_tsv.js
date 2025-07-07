function loadTSVTable(tsvPath, containerId) {
    console.log("Running loadTSVTable with:", tsvPath, containerId);
    fetch(tsvPath)
        .then(response => response.text())
        .then(tsv => {
            const rows = tsv.trim().split('\n').map(line => line.split('\t'));
            const table = document.createElement('table');
            table.classList.add('custom-table');

            const thead = document.createElement('thead');
            const headRow = document.createElement('tr');
            rows[0].forEach(cell => {
                const th = document.createElement('th');
                th.textContent = cell;
                headRow.appendChild(th);
            });
            thead.appendChild(headRow);
            table.appendChild(thead);

            const tbody = document.createElement('tbody');
            rows.slice(1).forEach(row => {
                const tr = document.createElement('tr');
                row.forEach(cell => {
                    const td = document.createElement('td');
                    td.textContent = cell;
                    tr.appendChild(td);
                });
                tbody.appendChild(tr);
            });
            table.appendChild(tbody);

            document.getElementById(containerId).appendChild(table);
        })
        .catch(error => {
            console.error("Error loading TSV:", error);
        });
}


function loadColorLegend(tsvPath, containerId) {
    fetch(tsvPath)
        .then(response => response.text())
        .then(tsv => {
            const lines = tsv.trim().split('\n');
            const headers = lines[0].split('\t');
            const rows = lines.slice(1).map(line => line.split('\t'));

            const container = document.getElementById(containerId);
            const legend = document.createElement('div');
            legend.classList.add('color-legend');

            rows.forEach(row => {
                const name = row[0];
                const r = parseFloat(row[2]) * 255;
                const g = parseFloat(row[3]) * 255;
                const b = parseFloat(row[4]) * 255;
                const color = `rgb(${r}, ${g}, ${b})`;

                const swatch = document.createElement('div');
                swatch.classList.add('color-swatch');
                swatch.style.backgroundColor = color;
                swatch.title = name;

                const label = document.createElement('span');
                label.textContent = name;

                const item = document.createElement('div');
                item.classList.add('color-item');
                item.appendChild(swatch);
                item.appendChild(label);

                legend.appendChild(item);
            });

            container.appendChild(legend);
        });
}