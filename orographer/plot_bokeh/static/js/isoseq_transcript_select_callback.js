export default function (args) {
const transcript_select = args.transcript_select;
const transcript_sources = args.transcript_sources;
const selectedTranscript = transcript_select.value || "ALL";
const selectedParts = String(selectedTranscript).split("::");
const selectedAnnotation = selectedParts.length == 2 ? selectedParts[0] : "";
const selectedTranscriptId = selectedParts.length == 2 ? selectedParts[1] : selectedTranscript;

function setSelectedIndices(source, indices) {
    source.selected.indices = indices;
    source.selected.change.emit();
}

function rowMatches(data, rowIndex) {
    const transcriptId = data.transcript_id ? data.transcript_id[rowIndex] : "";
    const annotationId = data.annotation_id ? data.annotation_id[rowIndex] : selectedAnnotation;
    if (String(transcriptId) != String(selectedTranscriptId)) {
        return false;
    }
    if (!selectedAnnotation) {
        return true;
    }
    return String(annotationId) == String(selectedAnnotation);
}

transcript_sources.forEach(function (source) {
    const data = source.data;
    const matches = [];
    if (selectedTranscript != "ALL") {
        (data.transcript_id || []).forEach(function (_transcriptId, rowIndex) {
            if (rowMatches(data, rowIndex)) {
                matches.push(rowIndex);
            }
        });
    }
    setSelectedIndices(source, matches);
    source.change.emit();
});
}
